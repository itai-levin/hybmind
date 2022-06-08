import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from joblib import Parallel, delayed

def reaction_ring_delta(rxn_smiles, retro=True):
    r, _, p = list(map(Chem.MolFromSmiles, rxn_smiles.split('>')))
    if retro:
        return rdMolDescriptors.CalcNumRings(p) - rdMolDescriptors.CalcNumRings(r)
    else:
        return rdMolDescriptors.CalcNumRings(r) - rdMolDescriptors.CalcNumRings(p)

def mol_num_chiral(mol):
    return sum([
        bool(atom.GetChiralTag())
        for atom in mol.GetAtoms()
    ])

def reaction_chiral_delta(rxn_smiles, retro=True):
    r, _, p = list(map(Chem.MolFromSmiles, rxn_smiles.split('>')))
    if retro:
        return mol_num_chiral(p) - mol_num_chiral(r)
    else:
        return mol_num_chiral(r) - mol_num_chiral(p)

def reaction_chiral_delta(rxn_smiles, retro=True):
    r, _, p = list(map(Chem.MolFromSmiles, rxn_smiles.split('>')))
    r.UpdatePropertyCache(strict=False)
    p.UpdatePropertyCache(strict=False)
    if retro:
        return rdMolDescriptors.CalcNumAtomStereoCenters(p) - rdMolDescriptors.CalcNumAtomStereoCenters(r)
    else:
        return rdMolDescriptors.CalcNumAtomStereoCenters(r) - rdMolDescriptors.CalcNumAtomStereoCenters(p)

def unmap(smiles, canonicalize=True):
    mol = Chem.MolFromSmiles(smiles)
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
    return Chem.MolToSmiles(mol, isomericSmiles=canonicalize)

def smiles_to_fingerprint(smi, length=2048, radius=2, useFeatures=False, useChirality=True):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        raise ValueError('Cannot parse {}'.format(smi))
    fp_bit = AllChem.GetMorganFingerprintAsBitVect(
        mol=mol, radius=radius, nBits = length, 
        useFeatures=useFeatures, useChirality=useChirality
    )
    return np.array(fp_bit)

def smiles_list_to_fingerprints(smiles_list, fp_length=2048, fp_radius=2, nproc=1):
    if nproc > 1:
        fps = Parallel(n_jobs=nproc, verbose=1)(
            delayed(smiles_to_fingerprint)(smi, length=fp_length, radius=fp_radius) for smi in smiles_list
        )
        fps = np.array(fps)
    else:
        fps = np.array([smiles_to_fingerprint(smi, length=fp_length, radius=fp_radius) for smi in smiles_list])
    return fps

def fix_spectators(reaction_smiles):
    reactants, spectators, products = reaction_smiles.split('>')
    reactants = Chem.MolFromSmiles(reactants)
    spectators = Chem.MolFromSmiles(spectators)
    products = Chem.MolFromSmiles(products)

    if not reactants or not spectators or not products:
        return reaction_smiles

    reactant_fragments = []
    spectator_fragments = []
    product_fragments = []

    # find reactants mislabeled as spectators first
    # these are spectator molecules with atom maps numbers that appear in the product
    # unset atom map number in spectators if it doesn't appear in products
    product_maps = set(a.GetAtomMapNum() for a in products.GetAtoms())
    for a in spectators.GetAtoms():
        if a.GetAtomMapNum() and a.GetAtomMapNum() not in product_maps:
            a.SetAtomMapNum(0)

    # identify the spectator molecules that should be reactants - move these to reactants
    for frag in Chem.MolToSmiles(spectators).split('.'):
        if not frag: continue # skip empty strings
        if any([bool(a.GetAtomMapNum()) for a in Chem.MolFromSmiles(frag).GetAtoms()]):
            reactant_fragments.append(frag)
        else:
            spectator_fragments.append(frag)

    # reactants need to be reconstructed with newly discovered molecules
    reactants = Chem.MolFromSmiles('.'.join([Chem.MolToSmiles(reactants)]+reactant_fragments))
    reactant_fragments = []

    reactant_maps = set(a.GetAtomMapNum() for a in reactants.GetAtoms())
    product_maps = set(a.GetAtomMapNum() for a in products.GetAtoms())

    # unset atom map number in products if it doesn't appear in reactants
    for a in products.GetAtoms():
        if a.GetAtomMapNum() and a.GetAtomMapNum() not in reactant_maps:
            a.SetAtomMapNum(0)

    # unset atom map number in reactants if it doesn't appear in products
    for a in reactants.GetAtoms():
        if a.GetAtomMapNum() and a.GetAtomMapNum() not in product_maps:
            a.SetAtomMapNum(0)

    # identify mislabeled reactants that do not contribute heavy atom to product - move these to spectators
    for frag in Chem.MolToSmiles(reactants).split('.'):
        if any([bool(a.GetAtomMapNum()) for a in Chem.MolFromSmiles(frag).GetAtoms()]):
            reactant_fragments.append(frag)
        else:
            spectator_fragments.append(frag)

    # identify mislabeled products that don't come from reactants - move these to spectators
    for frag in Chem.MolToSmiles(products).split('.'):
        if any([bool(a.GetAtomMapNum()) for a in Chem.MolFromSmiles(frag).GetAtoms()]):
            product_fragments.append(frag)
        else:
            spectator_fragments.append(frag)

    return '.'.join(reactant_fragments) + '>' + '.'.join(spectator_fragments) + '>' + '.'.join(product_fragments)

def precursors_from_template(mol, template):
    precursors = set()
    results = template.RunReactants([mol])
    for res in results:
        res_smiles = '.'.join(sorted([Chem.MolToSmiles(m, isomericSmiles=True) for m in res]))
        if Chem.MolFromSmiles(res_smiles):
            precursors.add(res_smiles)
    return list(precursors)

def precursors_from_templates(target_smiles, templates, nproc=1):
    mol = Chem.MolFromSmiles(target_smiles)
    precursor_set_list = Parallel(n_jobs=nproc, verbose=1)(
        delayed(precursors_from_template)(mol, template) for template in templates
    )
    return list(set().union(*precursor_set_list))

def reaction_from_smarts(smarts):
    return AllChem.ReactionFromSmarts(smarts)

def templates_from_smarts_list(smarts_list, nproc):
    templates = Parallel(n_jobs=nproc, verbose=1)(
        delayed(reaction_from_smarts)(smarts) for smarts in smarts_list
    )
    return templates

def bond_edit_stats(rsmi):
    r, s, p = rsmi.split('>')
    rmol = Chem.MolFromSmiles(r)    
    pmol = Chem.MolFromSmiles(p)
    
    pbonds = []
    for bond in pmol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a or b:
            pbonds.append(tuple(sorted([a, b])))
    
    rbonds = []
    for bond in rmol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a or b:
            rbonds.append(tuple(sorted([a, b])))
    
    r_changed = set(rbonds) - set(pbonds)
    p_changed = set(pbonds) - set(rbonds)
    atoms_changed = set()
    for ch in list(r_changed)+list(p_changed):
        atoms_changed.add(ch[0])
        atoms_changed.add(ch[1])
    atoms_changed -= set([0])
    return {
        'r_bonds': len(r_changed), 
        'p_bonds': len(p_changed), 
        'atoms': len(atoms_changed)
    }
