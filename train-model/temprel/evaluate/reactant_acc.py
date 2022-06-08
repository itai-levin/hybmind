import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRunText
from rdchiral.main import rdchiralRun
from joblib import Parallel, delayed

# int, Bool true_reactant_rank : For input product, reactant, k,
# and template-ranking vector, outputs the rank of the correct reactant and
# True if the reactant is in the top-k, else returns np.inf, False


def true_reactant_rank(prod_smiles, templates,
                       ranked_template_ids, true_reactant, k=1):

    # retrieve and unmap reactant
    react_mol = Chem.MolFromSmiles(true_reactant)
    for atom in react_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    true_reactant = Chem.MolToSmiles(react_mol)

    # rank the templates by highest likelihood
    ranked_template_ids = np.flip(np.argsort(ranked_template_ids))

    seen_reactants = set()
    n = 0

    # keep checking for true reactant until k unique reactants are found
    while len(seen_reactants) < k and n < len(ranked_template_ids):
        template_id = ranked_template_ids[n]
        pred_rxn = templates['reaction_smarts'][template_id]
        try:
            pred_reactants = rdchiralRunText(pred_rxn, prod_smiles)
        except:
            pred_reactants = []
        if true_reactant in pred_reactants:
            return n
        for reactant in pred_reactants:
            seen_reactants.add(reactant)
        n += 1

    return np.inf

# Float reactant_top_k_accuracy: For input dataframe with k, products, reactants,
# associated template-ranking vector, outputs the fraction of True/total entries
# when true_reactant_rank is mapped onto the dataframe.


def reactant_top_k_accuracy(ranks, k):
    # predict template rankings from input data
    length = len(ranks)
    results = (ranks <= (k-1))
    return np.sum(results) / length


def reactant_accuracy_by_popularity(model, test_smiles, pred_test_labels,
                                    test_labels, test_reactant_smiles, templates,
                                    train_labels, ks=[1, 10, 50, 100],
                                    njobs=2):
    class_counts = np.bincount(train_labels)  # number of training precedents
    cls_ind = {}
    for n in range(0, 11):
        cls_ind[n] = np.isin(test_labels, np.argwhere(
            class_counts == n).reshape(-1))
    cls_ind['med'] = np.isin(test_labels, np.argwhere(
        (class_counts > 10) & (class_counts <= 50)).reshape(-1))
    cls_ind['high'] = np.isin(
        test_labels, np.argwhere(class_counts > 50).reshape(-1))
    max_k = np.max(ks)
    ranks = Parallel(n_jobs=njobs)(delayed(true_reactant_rank)(
        test_smiles[i], templates, pred_test_labels[i], test_reactant_smiles[i], k=max_k)
        for i in range(len(pred_test_labels)))
    ranks = np.array(ranks)
    performance = {}
    for k in ks:
        performance[str(k)] = {}
        performance[str(k)]['all'] = reactant_top_k_accuracy(ranks, k)
        for n in range(0, 11):
            if len(pred_test_labels[cls_ind[n]]) == 0:
                continue
            print(n)
            performance[str(k)][n] = reactant_top_k_accuracy(
                ranks[cls_ind[n]], k)

        print('med')
        performance[str(k)]['med'] = reactant_top_k_accuracy(
            ranks[cls_ind['med']], k)
        print('high')
        performance[str(k)]['high'] = reactant_top_k_accuracy(
            ranks[cls_ind['high']], k)
    return performance
