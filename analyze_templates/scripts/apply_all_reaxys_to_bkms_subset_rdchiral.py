from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import re
from tqdm import tqdm
from rdchiral.main import rdchiralRun, rdchiralRunText
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
import datetime
import numpy as np

import os, contextlib

import warnings
warnings.filterwarnings("ignore")

import json

import argparse

def get_args():
    options = argparse.ArgumentParser()
    options.add_argument("--start", dest="start", type=int)
    options.add_argument("--end", dest="end", type=int)
    options.add_argument("--task-id", dest="task_id",type=str)
    options.add_argument("--reactions", dest="reactions",type=str)
    options.add_argument("--templates", dest="templates",type=str)
    options.add_argument("--prefix", dest="prefix",type=str,default=str(datetime.datetime.now()))
    args = options.parse_args()
    return args

def any_mapped_neighbors (atom):
    for neigh in atom.GetNeighbors():
        if neigh.GetAtomMapNum() != 0:
            return True
    return False

def neutralize_mol(mol, pattern):
    #pulled from http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
    """
    mol (rdkit.Chem.Mol): molecule to neutralize
    pattern (rdkit.Chem.Mol): pattern to match for neutralization - should be 
        defined as: pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        making it as a parameter to avoid re-instantiating 
    """
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
    return mol

def unmap(smarts):
    unmapped = re.sub('\:[0-9]+\]',']',smarts)
    return unmapped

def remove_degree(smarts):
    unmapped = re.sub('\&D[0-9]','',smarts)
    return unmapped

def remove_HCount(smarts):
    unmapped = re.sub('\&H[0-9]','',smarts)
    return unmapped

def show_reaction_list(list_of_reactions, smiles=True):
    if smiles:
        for rxn in list_of_reactions:
            try:
                print (rxn)
                plt.figure()
                plt.imshow(Chem.Draw.ReactionToImage(AllChem.ReactionFromSmarts(rxn, useSmiles=True)))
                plt.xticks([])
                plt.yticks([])
                plt.show()
            except Exception as e:
                print("Could not display reaction : ", rxn, e )
    else:
        for rxn in list_of_reactions:
            try:
                print (rxn)
                plt.figure()
                plt.imshow(Chem.Draw.ReactionToImage(AllChem.ReactionFromSmarts(rxn)))
                plt.xticks([])
                plt.yticks([])
                plt.show()
            except Exception as e:
                print("Could not display reaction : ", rxn, e )


def std_smarts_rxn (rxn_smarts):
    r, i, p = rxn_smarts.split('>')
    try:
        std_r = '.'.join(sorted(r.split('.')))
    except:
        std_r = r

    try:
        std_p = '.'.join(sorted(p.split('.')))
    except:
        std_p = p

    return '>'.join([std_r, i, std_p])

def flatten (ls):
    return [x for y in ls for x in y]

def rxns_recover_reactants(args, smarts=None, 
                            neut_pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")):
    idx, product, reactant = args
    for reaction_smarts in smarts:
        try:
            outcomes = [Chem.MolToSmiles(neutralize_mol(Chem.MolFromSmiles(x), neut_pattern),isomericSmiles=False) for x in rdchiralRun(reaction_smarts, product, combine_enantiomers=True)]
            for outcome in outcomes:
                if sorted(outcome.split(".")) == reactant:
                    return idx, True
        except (RuntimeError, ValueError, KeyError) as e:
            print (e)
            pass

    return idx, False

def get_rdchiral_reaction(smarts):
    try:
        return rdchiralReaction('('+smarts.replace('>>',')>>'))
    except:
        print ('_________FAILED___________')
        return rdchiralReaction('[C:1]>>[C:1]')

if __name__=='__main__':
    
    args = get_args()

    reactions = pd.read_json(args.reactions)
    all_templates = pd.read_json(args.templates)
    
    end = min(args.end, len(reactions))
    start = max(args.start, 0)
    
    reactions = reactions.iloc[start:end, :] # subset df

    print ('Getting product objects...')
    reactions['p_mols'] = [rdchiralReactants(Chem.MolToSmiles(Chem.MolFromSmiles(unmap(x)), isomericSmiles=False)) for x in reactions['products']]
    print ('Done getting product objects')

    print ('Getting reactant objects...')
    reactions['r_std_no_stereo_smiles'] = [sorted(Chem.MolToSmiles(Chem.MolFromSmiles(unmap(x)), isomericSmiles=False).split('.')) for x in reactions['reactants']]
    print ('Done getting reactant objects')
    print ('Getting reaction objects...', datetime.datetime.now())
    
    reaxys_smarts = []
    counter = 0
    smarts_list = all_templates[all_templates['template_set']=='reaxys']['reaction_smarts'].values

    for smarts in tqdm(smarts_list):
        reaxys_smarts.append(get_rdchiral_reaction(smarts))

    print ('Done getting reaction objects, {} reactions encoded,'.format(len(reaxys_smarts)))
    print ('Applying templates...', datetime.datetime.now())    


    recovered = []
    neut_pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    for idx, r, p in tqdm(zip(list(reactions.index), list(reactions['p_mols'].values), list(reactions['r_std_no_stereo_smiles'].values))):
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull): #hide output
                recovered.append(rxns_recover_reactants((idx, r, p), smarts=reaxys_smarts, neut_pattern=neut_pattern))


    print ('Done applying templates,', datetime.datetime.now())
    print (len(recovered))
    print (recovered[0])
    
    with open(args.prefix+'recovered_bkms_reactions_rdchiral_neutralized_{}.json'.format(str(args.task_id)), 'w') as f:
        json.dump(recovered, f)

