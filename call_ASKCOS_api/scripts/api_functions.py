import requests
import time
from pprint import pprint
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import json
import urllib3
urllib3.disable_warnings()


def get_all_children(tree, field = 'smiles', reactions=True):
"""
Paramaters:
tree: 
"""
    if not isinstance (tree, dict):
        print ("Input is not a dictionary, cannot process")
        return None
    reactions_list = []
    
    if reactions==True :
        query = 'is_reaction'
    else:
        query = 'is_chemical'
    try: 
        if tree[query]:
            reactions_list.append(tree[field])
            for child in tree['children']:
                reactions_list =   reactions_list + get_all_children(child, field, reactions)
    except:
        for child in tree['children']:
            reactions_list =  reactions_list + get_all_children(child, field, reactions)
    return reactions_list


