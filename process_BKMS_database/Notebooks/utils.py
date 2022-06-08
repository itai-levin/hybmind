import pandas as pd
import numpy as np
import re
import json
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse

from rdchiral.main import rdchiralRunText, rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants

#general helper functions

def standardize_smiles (smiles):
    mol = None
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
    if mol != None:
        return Chem.MolToSmiles(mol).replace('.[Na+]', '')
    else:
        return None


def mergeDict(dic1, dic2):
    """ 
    Merge dictionaries and keep values of common keys in list
    """
    dic = {}
    int_keys = set(dic1).intersection(set(dic2))
    dic3 = {key : list(set(dic1[key]).union(dic2[key])) for key in int_keys}
    dic.update(dic1)
    dic.update(dic2)
    dic.update(dic3)
    return dic

def intersection (l1, l2):
    """
    Returns the intersection of two collections, l1 and l2
    """
    return [i for i in l1 if i in l2]

def get_num_products (reaction_smiles):
    
    if reaction_smiles:
        return len(reaction_smiles.split('>>')[1].split('.'))
    else:
        return reaction_smiles
    
def neutralize_atoms(smiles):
    """
    from http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    
    input:
        smiles (str): molecule SMILES string
    
    output:
        SMILES string of input molecule 
    """
    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return Chem.MolToSmiles(mol)

def neutralize_reaction (reaction_smiles):
    try:
        return '>>'.join([neutralize_atoms(smiles) for smiles in reaction_smiles.split('>>')])
    except:
        print ("Couldn't neutralize,", reaction_smiles)
        return reaction_smiles

def show_reaction(reaction_smiles, useSmiles=True, subImgSize=(1200,1200)):
    plt.figure(figsize=(20,20))
    plt.imshow(Chem.Draw.ReactionToImage(AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True), subImgSize=(1200,1200)))
    plt.xticks([])
    plt.yticks([])
    plt.show()

#functions to map from BKMS database names to SMILES

def list_to_smiles (ls, dic, standardize = False) :
    """
    Converts a list of common chemicals to a list of SMILES using a dictionary 
    to translate between the two
    
    parameters: 
        ls (list[str]): list of chemical names 
        dic (dict): dictionary of the form {chemical name:SMILES string}
    
    output: 
        string of reaction SMILES
    
    note:
        if any chemical unrecognized, returns None
    """
    smiles_ls = []
    for i in ls:
        s = re.sub("^[Aa]n? |^n ", "", i).strip() #remove articles (a, an..., and ambiguous coefficients)
        #handle stoichiometric coeffs
        try:
            coeff = re.match("^[0-9]+ ", i)
            if coeff != None:
                coeff = coeff.group().strip()
                if coeff.isdigit():
                    s = re.sub("[0-9]+ ", "", s)
                    smi = dic[s]
                    if standardize:
                        smi = standardize_smiles(smi)
                    smiles_ls.extend([smi] * int(coeff))

            else:
                smi = dic[s]
                if standardize:
                    smi = standardize_smiles(smi)
                smiles_ls.append(smi)
        except (AttributeError, KeyError):
            return None
    return '.'.join(smiles_ls)


def conv_to_smiles (string, dic_ls, cofs, lower = True, debugging=False) :
    """
    parameters: 
        string (str): reaction string in the form "chem_1 + chem_2 + ... = chem_3 + chem_4 + ..."
        dic_ls (list[dict]): dictionary with chemical name strings as keys and SMILES strings as values
        cofs (Iterable): set of chemical names of "cofactors"

    output: 
        reaction as SMILES string
    """
    if lower:
        string = string.lower()
    if not isinstance(dic_ls, list):
        dic_ls = [dic_ls]
    split = re.split(" = | <=> ", string)
    if len(split) == 1: return dic_ls[0][split[0]]

    simp = re.compile("^[Aa]n? |^[0-9]+ |^n ")
    
    reactants = re.split(" \+ ", split[0])
    products = re.split(" \+ ", split[1])
    
    reactant_keys = [simp.sub("", reactant).lower() for reactant in reactants]
    product_keys = [simp.sub("", product).lower() for product in products]
    
    r = []
    p = []
    try: #if chemicals to be ignored are specified
        r_rm = list(cofs["ignore"]) #list of reactants which should be removed based on the products
        p_rm = list(cofs["ignore"]) #list of products which should be removed based on the reactants
    except KeyError:
        r_rm = [] #list of reactants which should be removed based on the products
        p_rm = [] #list of products which should be removed based on the reactants
    
    #create list of chemical names to be removed from the reaction string
    for key in reactant_keys:
        if key in cofs.keys() : p_rm.extend(cofs[key])
    for key in product_keys:
        if key in cofs.keys() : r_rm.extend(cofs[key])
    
    for key, reactant in zip(reactant_keys, reactants):
        if (key not in cofs.keys() and key not in r_rm) or\
            (key in cofs.keys() and\
             intersection(cofs[key], product_keys) == [] ) :
            r.append(reactant)
    for key, product in zip(product_keys, products):
        if (key not in cofs.keys() and key not in p_rm) or\
            (key in cofs.keys() and\
             intersection(cofs[key], reactant_keys) == [] ) :
            p.append(product)
    if len(r) == 0 or len(p) == 0: #if removing cofactors elimnates all products or reactants, restore them
        r = reactants
        p = products

    #to increase database consistency, keeping the dicts separate and try each one individually
    for dic in dic_ls:
        r_smiles = list_to_smiles(r, dic)
        p_smiles = list_to_smiles(p, dic)
        if r_smiles != None and p_smiles != None:
            break
    
    if r_smiles == None or p_smiles == None:
        if debugging:
            return str(r_smiles) + ">>" + str(p_smiles)
        else:
            return None
        
    else:
        return r_smiles+">>"+p_smiles

    
#functions to automatically identify co-factors

def get_prods_reactants (rxn_str_list, f_or_r):
    """
    Returns a set of string for the reactants, products, and reactions using chemical names
    
    Parameters:
    rxn_str_list list[str]: reactions with chemical names
    f_or_r (char): 'f' or 'r' whether to look at rxns in the fwd or reverse direction
    
    Outputs:
    rxtnts (set[string]): of reactant chemical names
    prods (set[string]): product chemical names
    rxns (set[string]): reactions
    """
    prods = set([])
    rxtants = set([])
    rxns = set([])
    simp = re.compile("^[Aa]n? |^[0-9]+ |^n ")
    for rxn_str in rxn_str_list:
        rxns.add(rxn_str.lower())
        split1 = re.split(" = | <=> ", rxn_str.lower())

        if f_or_r == "f":
            reactants = re.split(" \+ ", split1[0])
            products = re.split(" \+ ", split1[1])
        elif f_or_r == "r":
            reactants = re.split(" \+ ", split1[1])
            products = re.split(" \+ ", split1[0])
        else:
            print ("forward or backward reactions not specified")
            raise
        
        for x in reactants:
            s = simp.sub( "", x) #remove articles/ stoichiometric coefficients
            rxtants.add(s)
        for x in products:
            s = simp.sub("", x)
            prods.add(s)
    return rxtants, prods, rxns



def make_coOccurence_tab (rxn_str_list, f_or_r):
    """
    
    Parameters:
        rxn_str_list list[string]: reactions of form 'Chem_1 [+ Chem_2 + ... ]= Chem 3 ...'
        f_or_r (string): 'f' or 'r' indicates whether to look at rxns in the fwd or reverse direction
    
    Output:
        numpy array (size num_reactants*num_products) where each entry corresponds to the number of 
        reactions in which a  pair of molecules co-occur as a reactant-product pair
    
    notes: raises an exception if any value other than 'f' or 'r' is passes as f_or_r
    """
    all_reactants, all_products, all_reactions = get_prods_reactants(rxn_str_list, f_or_r)
    all_reactants = list(all_reactants)
    all_products = list(all_products)
    all_reactions = list(all_reactions)
    reactant_to_ind = {}
    product_to_ind = {}
    tot_occurences = np.zeros(len(all_reactants))
    simp = re.compile("^[Aa]n? |^[0-9]+ |^n ")
    for i in range(len(all_reactants)):
        reactant_to_ind[all_reactants[i]] = i
    for i in range(len(all_products)):
        product_to_ind[all_products[i]] = i
    tab = np.zeros((len(all_reactants),len(all_products)))
    
    for rxn_str in all_reactions:
        split1 = re.split(" = | <=> ", rxn_str)
        if f_or_r == "f":
            reactants = re.split(" \+ ", split1[0])
            products = re.split(" \+ ", split1[1])
        elif f_or_r == "r":
            reactants = re.split(" \+ ", split1[1])
            products = re.split(" \+ ", split1[0])
        else:
            print ("forward or backward reactions not specified")
            raise
        for reactant in reactants:
            r = simp.sub("", reactant).lower() #remove articles/ stoichiometric coefficients
            r_ind = reactant_to_ind[r]
            tot_occurences[r_ind] += 1
            for prod in products:
                p = simp.sub("", prod).lower() #remove articles/ stoichiometric coefficients
                p_ind = product_to_ind[p]
                tab[r_ind, p_ind] += 1

            
    return tab, tot_occurences, all_reactants, all_products, reactant_to_ind, product_to_ind


def make_cofactor_dict(rxn_str_list, occ_cutoff, frac_cutoff, f_or_r):
    """
    Returns a dictionary of 'cofactors,' defined by their appearance in over [occ_cutoff]
    reactions and coappearance with another molecule over [frac_cutoff] of the time
    
    Parameters:
        rxn_str_list (list[string]): reactions represented with common chemical names
        occ_cutoff (int): minimum number of appearances
        frac_cutoff (float): minimum fraction of coappearances
        f_or_r (string): 'f' or 'r' whether to look at rxns in the fwd or reverse direction
    
    output: 
        dict of form :{chemical name: coappearing chemical name}
    """
    
    tab, tot_occurences, reactants, products, r_dic, p_dic = make_coOccurence_tab (rxn_str_list, f_or_r)
    r = set()
    p = {}
    rel_tab = tab/tot_occurences[:,None]
    #indices for reactants that appear more than 30 times
    for i in np.arange(len(tot_occurences))[tot_occurences > occ_cutoff]: 
        chem = reactants[i]
        if np.max(rel_tab [r_dic[chem],:]) > frac_cutoff:
            r.add(reactants[i])
            l = []
            for j in np.arange(len(rel_tab [r_dic[chem],:])) [np.argsort(rel_tab [r_dic[chem],:])][np.sort(rel_tab [r_dic[chem],:])>frac_cutoff]:
                l.append(products[j])
            p[reactants[i]] = l
    return p
                   

# functions used to process mapped SMILES

def remove_molecule(query_mol, smiles):
    """
    Parameters:
        query_mol (rdkit.Mol): molecule to remove from a molecule SMILES
        smiles (str): molecule SMILES from which to remove molecule
    
    output:
        SMILES with query_mol removed if it matched one of the molecules 
    """
    new_mols_list = []
    for smile in smiles.split('.'):
        mol = Chem.MolFromSmiles(smile)
        if mol:
            if not mol.HasSubstructMatch(query_mol) or not query_mol.HasSubstructMatch(mol):
                new_mols_list.append(mol)
        else:
            return smiles
    return '.'.join([Chem.MolToSmiles(mol) for mol in new_mols_list])

def remove_molecule_from_reactants_and_products (query_mol, reaction_smiles, products_only=False, reactants_only=False):
    """
    Parameters:
        query_mol (rdkit.Mol): molecule to remove from a molecule SMILES
        reaction_smiles (str): reaction SMILES from which to remove moleule
        products_only (bool): if True, remove the molecule only from the products side
        reactants_only (bool):if True, remove the molecule only from the reactants side
    
    Output
    """
    assert not (products_only and reactants_only)
    if '>' not in reaction_smiles:
        return None
    reactants, intermediate, products = reaction_smiles.split('>')
    if products_only:
        return '>'.join([reactants, intermediate, remove_molecule(query_mol, products)])
    elif reactants_only:
        return '>'.join([remove_molecule(query_mol, reactants), intermediate, products])
    return '>'.join([remove_molecule(query_mol, reactants),intermediate, remove_molecule(query_mol, products)])

def unmap_unmatched_products(smiles):
    """
    Remove the atom map number for atoms on the products if they are not matched to 
    atoms in the reactants
    """
    reactants, i, products = smiles.split('>')
    prod = Chem.MolFromSmiles(products)
    react_map_nums = [at.GetAtomMapNum() for at in Chem.MolFromSmiles(reactants).GetAtoms()]
    [at.SetAtomMapNum(0) for at in prod.GetAtoms() if at.GetAtomMapNum() not in react_map_nums]
    return '>'.join([reactants, i, Chem.MolToSmiles(prod)])

def unmap_unmatched_reactants(smiles):
    """
    Remove the atom map number for atoms on the reactant if they are not matched to 
    atoms in the products
    """
    reactants, i, products = smiles.split('>')
    react = Chem.MolFromSmiles(reactants)
    prod_map_nums = [at.GetAtomMapNum() for at in Chem.MolFromSmiles(products).GetAtoms()]
    [at.SetAtomMapNum(0) for at in react.GetAtoms() if at.GetAtomMapNum() not in prod_map_nums]
    return '>'.join([Chem.MolToSmiles(react), i, products])

def flip_reaction (reaction):
    try:
        if '<=>' in reaction: # BKMS reaction string
            reactants, products = reaction.split(' <=> ')
            return ' <=> '.join([products, reactants])
        else: # smiles
            reactants, i, products = reaction.split('>')
            return '>'.join([products, i, reactants])
    except ValueError:
        return None

def reverse_reactions_in_df (df):
    df['smiles'] = df['smiles'].map(lambda x: flip_reaction(x))
    df['mapped_smiles'] = df['mapped_smiles'].map(lambda x: flip_reaction(x))


def count_second_product(reaction_ls):
    products = [x.split('>')[2] for x in reaction_ls if len(x.split('>')[2].split('.'))>1]
    freq_dict = {}
    for p in products:
        mols = p.split('.')
        for mol in mols:
            if mol in freq_dict.keys():
                freq_dict[mol] += 1
            else:
                freq_dict[mol] = 1
    return freq_dict

def num_mapped (mol):
    mol = Chem.MolFromSmiles(mol)
    mapped = 0
    total = 0
    for at in mol.GetAtoms():
        total += 1
        if at.GetAtomMapNum() != 0:
            mapped += 1
    return mapped, total


def remove_poorly_mapped_reactants(smiles, thresh = 0.15, min_mapped = 3):
    """
    smiles: string 
        reaction SMILES
    thresh: float 
        Threshold fraction of atoms required to keep a molecule
    min_mapped: int
        Minimum number of atoms required to keep a molucule
    """
    reactants, i, products = smiles.split('>')
    reactants = reactants.split('.')
    if len(reactants) > 1:
        mapped, total = zip(*[num_mapped(mol) for mol in reactants])
        frac_mapped = {smile:(m/t) for smile, m, t in zip(reactants, mapped, total)}
        abs_mapped = {smile:m for smile,m in zip(reactants, mapped)}
        new_reactants = []
        for smile in frac_mapped:
            if frac_mapped[smile] > thresh and abs_mapped[smile] >= min_mapped:
                new_reactants.append(smile)
        return '>'.join(['.'.join(new_reactants), i, products])
    else:
        return smiles

def has_poorly_mapped_reactants(smiles, thresh = 0.15, min_mapped = 3):
    """
    smiles: string 
        reaction SMILES
    thresh: float 
        Threshold fraction of atoms required to keep a molecule
    min_mapped: int
        Minimum number of atoms required to keep a molucule
    """
    reactants, i, products = smiles.split('>')
    reactants = reactants.split('.')
    if len(reactants) > 1:
        mapped, total = zip(*[num_mapped(mol) for mol in reactants])
        frac_mapped = {smile:(m/t) for smile, m, t in zip(reactants, mapped, total)}
        abs_mapped = {smile:m for smile,m in zip(reactants, mapped)}
        new_reactants = []
        for smile in frac_mapped:
            if frac_mapped[smile] < thresh or abs_mapped[smile] < min_mapped:
                return True
        return False
    else:
        return False

def same_products_and_reactants_smiles (smiles):
    return smiles.split('>')[0] == smiles.split('>')[2]

def same_products_and_reactants_names (reaction):
    return reaction.split('=')[0] == reaction.split('=')[1]
