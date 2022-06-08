from joblib import Parallel, delayed
from ..rdkit import reaction_ring_delta, reaction_chiral_delta

def calculate_attributes(df, nproc=1):
    df['ring_delta'] = Parallel(n_jobs=nproc, verbose=1)(
        delayed(reaction_ring_delta)(rxnsmiles) for rxnsmiles in df['reaction_smiles']
    )
    df['chiral_delta'] = Parallel(n_jobs=nproc, verbose=1)(
        delayed(reaction_chiral_delta)(rxnsmiles) for rxnsmiles in df['reaction_smiles']
    )
