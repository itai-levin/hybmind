import os
import re
import json
import hashlib
import numpy as np
import pandas as pd
import tqdm
from joblib import Parallel, delayed
from multiprocessing import Pool, TimeoutError
from .filter import filter_by_bond_edits
from .attributes import calculate_attributes
from .validate import validate_template
from ..rdkit import smiles_to_fingerprint, fix_spectators
from ..rdkit import unmap as rdkit_unmap

import sys
sys.path.append('/Users/Itai/Box Sync/Grad/research/rdchiral')
import rdchiral
from rdchiral.template_extractor import extract_from_reaction
print(rdchiral.template_extractor.__file__)

def create_hash(pd_row):
    return hashlib.md5(pd_row.to_json().encode()).hexdigest()


def unmap(smarts):
    return re.sub(r":[0-9]+]", "]", smarts)


def count_products(smiles):
    return 1 + smiles.count(".")


def extract_template(reaction):
    try:
        return extract_from_reaction(reaction)
    except KeyboardInterrupt:
        print("Interrupted")
        raise KeyboardInterrupt
    except Exception as e:
        return {"reaction_id": reaction["_id"], "error": str(e)}


def parallel_extract(df, nproc, timeout=3):
    pool = Pool(processes=nproc)
    async_results = [
        pool.apply_async(extract_template, (r,)) for r in df.to_dict(orient="records")
    ]
    templates = []
    for res in tqdm.tqdm(async_results):
        try:
            templates.append(res.get(timeout))
        except TimeoutError:
            templates.append(None)
    return templates


def templates_from_reactions(
    df,
    output_prefix=None,
    nproc=8,
    timeout=3,
    filter_fn=filter_by_bond_edits,
    filter_kwargs={},
):

    assert output_prefix is not None, "Please specify a directory for the output files"
    print ('INITIAL', df[df['EC_Number'].map(lambda x : '4.3.1.25' in str(x))])
    if not os.path.exists(output_prefix):
        os.makedirs(output_prefix)

    out_json = output_prefix + "templates.df.json.gz"

    assert type(df) is pd.DataFrame
    assert "reaction_smiles" in df.columns
    df["reaction_smiles"] = Parallel(n_jobs=nproc, verbose=1)(
        delayed(fix_spectators)(rsmi) for rsmi in df["reaction_smiles"]
    )
    rxn_split = df["reaction_smiles"].str.split(">", expand=True)
    df[["reactants", "spectators", "products"]] = rxn_split.rename(
        columns={0: "reactants", 1: "spectator", 2: "products"}
    )
    if "_id" not in df.columns:
        df["_id"] = df.apply(create_hash, axis=1)
    num_products = df["products"].apply(count_products)
    df = df[num_products == 1]
    templates = parallel_extract(df, nproc, timeout)  # extract templates
    templates = pd.DataFrame(filter(lambda x: x, templates)) #filter out False templates 
    keep_cols = [
        "dimer_only",
        "intra_only",
        "necessary_reagent",
        "reaction_id",
        "reaction_smarts",
    ]
    df = df.merge(templates[keep_cols], left_on="_id", right_on="reaction_id")
    df = df.dropna(subset=["reaction_smarts"]) #remove NaN templates
    print("DF length", len(df))
    valid = Parallel(n_jobs=nproc, verbose=1)(
        delayed(validate_template)(template)
        for template in df.to_dict(orient="records")
    ) #validate
    df = df[valid]
    print("valid length", len(df))
    df = filter_fn(df, nproc=nproc, **filter_kwargs) #irrelevant filter
    print("Length after filtering", len(df))
    #count 
    df["unmapped_template"] = df["reaction_smarts"].apply(unmap)
    unique_templates = pd.DataFrame(
        df["unmapped_template"].unique(), columns=["unmapped_template"]
    ).reset_index()
    df = df.merge(unique_templates, on="unmapped_template")
    print("Length after merging", len(df))
    df = df.drop_duplicates(subset=["products", "index"])
    print("Length after dropping duplicates", len(df))
    df["count"] = df.groupby("index")["index"].transform("count")
    if out_json[-2] == "gz":
        df.to_json(out_json, compression="gzip")
    else:
        df.to_json(out_json)
    return df


def create_split(df, method, class_col="index"):
    df = df.sample(frac=1).reset_index(drop=True)
    if method == "random":
        i80 = int(len(df) * 0.8)
        i90 = int(len(df) * 0.9)
        df.loc[:i80, "split"] = "train"
        df.loc[i80:i90, "split"] = "valid"
        df.loc[i90:, "split"] = "test"
    elif method == "stratified":
        labels = df[class_col].values
        class_counts = np.bincount(df["index"].values)
        classes_with_1_example = np.argwhere(class_counts == 1).reshape(-1)
        np.random.shuffle(classes_with_1_example)
        classes_with_2_example = np.argwhere(class_counts == 2).reshape(-1)
        classes_with_few_example = np.argwhere(
            (class_counts > 2) & (class_counts <= 10)
        ).reshape(-1)
        classes_with_many_example = np.argwhere(class_counts > 10).reshape(-1)
        train_ind = []
        val_ind = []
        test_ind = []
        ind = np.argwhere(np.isin(labels, classes_with_1_example)).reshape(-1)
        i80 = int(len(ind) * 0.8)
        i90 = int(len(ind) * 0.9)
        train_ind.extend(ind[:i80])
        val_ind.extend(ind[i80:i90])
        test_ind.extend(ind[i90:])
        for cls in classes_with_2_example:
            ind = np.argwhere(np.isin(labels, cls)).reshape(-1)
            train_ind.append(ind[0])
            test_ind.append(ind[1])
        for cls in classes_with_few_example:
            ind = np.argwhere(np.isin(labels, cls)).reshape(-1)
            np.random.shuffle(ind)
            test_ind.append(ind[0])
            val_ind.append(ind[1])
            train_ind.extend(ind[2:])
        for cls in classes_with_many_example:
            ind = np.argwhere(np.isin(labels, cls)).reshape(-1)
            np.random.shuffle(ind)
            i80 = int(len(ind) * 0.8)
            i90 = int(len(ind) * 0.9)
            train_ind.extend(ind[:i80])
            val_ind.extend(ind[i80:i90])
            test_ind.extend(ind[i90:])
        df.loc[train_ind, "split"] = "train"
        df.loc[val_ind, "split"] = "valid"
        df.loc[test_ind, "split"] = "test"
    df = df.sample(frac=1)
    return df


def process_for_training(
    templates_df,
    output_prefix,
    split_col="split",
    calc_split=None,
    smiles_col="products",
    class_col="index",
    react_col="reactants",
):
    if calc_split in ["random", "stratified"]:
        templates_df = create_split(templates_df, calc_split, class_col)

    if split_col not in templates_df.columns:
        raise ValueError('split column "{}" not in DataFrame.')

    if smiles_col not in templates_df.columns:
        ValueError('smiles column "{}" not in DataFrame.')

    if class_col not in templates_df.columns:
        ValueError('class column "{}" not in DataFrame.')

    if output_prefix is None:
        output_prefix = ""

    train_df = templates_df[templates_df[split_col] == "train"]
    valid_df = templates_df[templates_df[split_col] == "valid"]
    test_df = templates_df[templates_df[split_col] == "test"]

    train_smiles = train_df[smiles_col].values.astype(str)
    valid_smiles = valid_df[smiles_col].values.astype(str)
    test_smiles = test_df[smiles_col].values.astype(str)

    train_classes = train_df[class_col].values.astype(int)
    valid_classes = valid_df[class_col].values.astype(int)
    test_classes = test_df[class_col].values.astype(int)

    train_reacts = train_df[react_col].values.astype(str)
    valid_reacts = valid_df[react_col].values.astype(str)
    test_reacts = test_df[react_col].values.astype(str)

    np.save(output_prefix + "train.input.smiles.npy", train_smiles)
    np.save(output_prefix + "valid.input.smiles.npy", valid_smiles)
    np.save(output_prefix + "test.input.smiles.npy", test_smiles)

    np.save(output_prefix + "train.labels.classes.npy", train_classes)
    np.save(output_prefix + "valid.labels.classes.npy", valid_classes)
    np.save(output_prefix + "test.labels.classes.npy", test_classes)

    np.save(output_prefix + "train.reactants.smiles.npy", train_reacts)
    np.save(output_prefix + "valid.reactants.smiles.npy", valid_reacts)
    np.save(output_prefix + "test.reactants.smiles.npy", test_reacts)


def process_for_askcos(templates_df, template_set_name, output_prefix, nproc=8):
    templates_df["reaction_id"] = range(len(templates_df))
    template_references = templates_df.groupby("index")["reaction_id"].apply(list)
    templates_df["references"] = templates_df["index"].map(template_references)
    templates_df["template_set"] = template_set_name
    template_columns_to_keep = [
        "index",
        "reaction_smarts",
        "necessary_reagent",
        "intra_only",
        "dimer_only",
        "count",
        "template_set",
        "references",
    ]
    templates_to_save = templates_df.sort_values("index").drop_duplicates("index")[
        template_columns_to_keep
    ]
    reactants_concat = templates_df["reactants"].values
    products = templates_df["products"].values
    reactants = []
    for react in reactants_concat:
        reactants.extend(react.split("."))
    reactants = Parallel(n_jobs=nproc, verbose=1)(
        delayed(rdkit_unmap)(smiles) for smiles in reactants
    )
    products = Parallel(n_jobs=nproc, verbose=1)(
        delayed(rdkit_unmap)(smiles) for smiles in products
    )
    react_df = pd.DataFrame(reactants, columns=["as_reactant"])
    prod_df = pd.DataFrame(products, columns=["as_product"])
    react_count = react_df.groupby("as_reactant")["as_reactant"].count()
    prod_count = prod_df.groupby("as_product")["as_product"].count()
    react_count = react_count.to_frame()
    react_count.index.name = "smiles"
    prod_count = prod_count.to_frame()
    prod_count.index.name = "smiles"
    historian_df = (
        react_count.merge(prod_count, how="outer", left_index=True, right_index=True)
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    historian_df["template_set"] = template_set_name

    calculate_attributes(templates_df, nproc=nproc)
    templates_to_save["ring_delta"] = templates_to_save["index"].map(
        templates_df.groupby("index")["ring_delta"].mean()
    )
    templates_to_save["chiral_delta"] = templates_to_save["index"].map(
        templates_df.groupby("index")["chiral_delta"].mean()
    )
    templates_to_save["attributes"] = templates_to_save[
        ["ring_delta", "chiral_delta"]
    ].to_dict(orient="records")
    templates_to_save = templates_to_save.drop(["ring_delta", "chiral_delta"], axis=1)

    templates_to_save["_id"] = templates_to_save.apply(create_hash, axis=1)
    templates_to_save.to_json(
        output_prefix + "retro.templates.{}.json.gz".format(template_set_name),
        orient="records",
        compression="gzip",
    )

    reactions_df = templates_df.drop(
        columns=template_columns_to_keep + ["unmapped_template"]
    )
    reactions_df["template_set"] = template_set_name
    reactions_df.to_json(
        output_prefix + "reactions.{}.json.gz".format(template_set_name),
        orient="records",
        compression="gzip",
    )

    historian_df.to_json(
        output_prefix + "historian.{}.json.gz".format(template_set_name),
        orient="records",
        compression="gzip",
    )
