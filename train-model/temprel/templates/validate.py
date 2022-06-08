import sys
sys.path.append('/Users/Itai/Box Sync/Grad/research/rdchiral')
import rdchiral
from rdchiral.main import rdchiralRunText
from rdkit import Chem
from rdkit.Chem import AllChem
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.main import rdchiralRunText
from rdchiral.main import rdchiralRun
import re

print(rdchiral.__file__)

def union(list_of_lists):
    u = set()
    for ls in list_of_lists:
        u = u.union(set(ls))
    return u


def intersection(list_of_lists):
    i = set(list_of_lists[0])
    for ls in list_of_lists[1:]:
        i = i.intersection(set(ls))
    return i

def smiles_processing(smiles):
    """
    apply processing function to smiles
    """
    smiles = re.sub('[/\\\]C=N[/\\\]', 'C=N', smiles)
    return smiles

def simplify_stereo(template):
    smarts = template["reaction_smarts"]
    p_smarts, _, r_smarts = smarts.split(">")
    
    reactants = Chem.MolFromSmiles(smiles_processing(template["reactants"]))
    products = Chem.MolFromSmiles(smiles_processing(template["products"]))
    

    r_smarts = AllChem.MolFromSmarts(r_smarts)
    p_smarts = AllChem.MolFromSmarts(p_smarts)

    # Find all atoms matched by the template to preserve stereo-information
    react_ctr_idxs = union(reactants.GetSubstructMatches(r_smarts, useChirality=True))
    prod_ctr_idxs = union(products.GetSubstructMatches(p_smarts, useChirality=True))

    react_ctr_maps = [
        reactants.GetAtomWithIdx(i).GetAtomMapNum() for i in react_ctr_idxs
    ]
    prod_ctr_maps = [products.GetAtomWithIdx(i).GetAtomMapNum() for i in prod_ctr_idxs]

    react_maps = intersection([react_ctr_maps, prod_ctr_maps])

    react_idxs = {
        at.GetAtomMapNum(): at.GetIdx()
        for at in reactants.GetAtoms()
        if at.GetAtomMapNum() not in react_maps
    }
    prod_idxs = {
        at.GetAtomMapNum(): at.GetIdx()
        for at in products.GetAtoms()
        if at.GetAtomMapNum() not in react_maps
    }

    mapnums = intersection((prod_idxs.keys(), react_idxs.keys()))

    for mapnum in mapnums:
        
        #remove chiral information from atoms not in reaction center
        distant_react_atom = reactants.GetAtomWithIdx(react_idxs[mapnum])
        distant_prod_atom = products.GetAtomWithIdx(prod_idxs[mapnum])
        distant_react_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
        distant_prod_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
        
        #remove E/Z information from bonds not in reaction center
        for bond in distant_react_atom.GetBonds():
            if bond.GetBeginAtom().GetAtomMapNum() in mapnums and bond.GetEndAtom().GetAtomMapNum() in mapnums:
                bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)
                bond.SetBondDir(Chem.rdchem.BondDir.NONE)
                
        for bond in distant_prod_atom.GetBonds():
            if bond.GetBeginAtom().GetAtomMapNum() in mapnums and bond.GetEndAtom().GetAtomMapNum() in mapnums:
                bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)
                bond.SetBondDir(Chem.rdchem.BondDir.NONE)


    return reactants, Chem.MolToSmiles(products)


def validate_template(template):

    reaction_smarts = template["reaction_smarts"]

    reactants, products = simplify_stereo(template)
    
    reaction_smarts = "(" + reaction_smarts.replace(">>", ")>>")  # make unimolecular
    try:
        outcomes = [Chem.MolToSmiles(Chem.MolFromSmiles(x),isomericSmiles=True) for x in rdchiralRunText(reaction_smarts, products, combine_enantiomers=False)]
    except (RuntimeError, ValueError, KeyError) as e:
        return False

    # Clear reactants map to get canonical SMILES
    for a in reactants.GetAtoms():
        a.SetAtomMapNum(0)

    #this seems excessive, but in practice seems to make a difference
    reactants_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(reactants, isomericSmiles=True)))

    for outcome in outcomes:
        if sorted(outcome.split(".")) == sorted(reactants_canonical.split(".")):
            return True
    return False


def mapping_validity(reaction_smiles):
    reactants, spectators, products = reaction_smiles.split(">")
    reactant_mol = Chem.MolFromSmiles(reactants)
    product_mol = Chem.MolFromSmiles(products)

    reactant_mapping = {}
    for atom in reactant_mol.GetAtoms():
        map_number = atom.GetAtomMapNum()
        if not map_number:
            continue
        if map_number in reactant_mapping:
            return "DUPLICATE_REACTANT_MAPPING"
        reactant_mapping[map_number] = atom.GetIdx()

    product_mapping = {}
    for atom in product_mol.GetAtoms():
        map_number = atom.GetAtomMapNum()
        if not map_number:
            continue
        if map_number in product_mapping:
            return "DUPLICATE_PRODUCT_MAPPING"
        product_mapping[map_number] = atom.GetIdx()

    if len(reactant_mapping) < len(product_mapping):
        return "UNMAPPED_REACTANT_ATOM(S)"

    for map_number in product_mapping.keys():
        if map_number not in reactant_mapping:
            return "UNMAPPED_PRODUCT_ATOM"
        reactant_atom = reactant_mol.GetAtomWithIdx(reactant_mapping[map_number])
        product_atom = product_mol.GetAtomWithIdx(product_mapping[map_number])

        if reactant_atom.GetSymbol() != product_atom.GetSymbol():
            return "ALCHEMY"

    return "VALID"
