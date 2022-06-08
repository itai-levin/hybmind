import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from ..rdkit import bond_edit_stats

def filter_by_bond_edits(df, max_rbonds=5, max_pbonds=3, max_atoms=10, nproc=1, return_filtered=False):
    """Filter templates by bond edit statistics.

    Args:
        max_rbonds (int): Maximum number of reactant bonds that can change (others get removed)
        max_pbonds (int): Maximum number of product bonds that can change
        max_atoms (int): Maximum number of atoms that can change connectivity
    
    Returns:
        pd.DataFrame: data frame of templates that survived the filter
        pd.DataFrame: if `return_filtered==True`, return a separate dataframe of the templates that were removed
    """
    bond_stats = Parallel(n_jobs=nproc, verbose=1)(
        delayed(bond_edit_stats)(rsmi) for rsmi in df['reaction_smiles']
    )

    bond_stats = pd.DataFrame(bond_stats)
    
    bool_mask = (
        (bond_stats['r_bonds'] <= max_rbonds) & 
        (bond_stats['p_bonds'] <= max_pbonds) & 
        (bond_stats['atoms'] <= max_atoms)
    )

    filtered = df.iloc[~bool_mask.values]

    df = df.iloc[bool_mask.values]

    if return_filtered:
        return df, filtered
    
    return df