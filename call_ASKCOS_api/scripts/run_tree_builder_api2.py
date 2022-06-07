import requests
import time
from pprint import pprint
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import json
import urllib3
from joblib import Parallel, delayed

urllib3.disable_warnings()

def save_results (results, smiles, path = None, write_type='w'):
    if path == None:
        path = 'output_trees/{}_paths.json'.format(smiles)
    with open(path, write_type) as outfile:
        json.dump(results, outfile)
        outfile.write('\n')

def post_and_get (HOST, params, sleep_time = 10, timeout = 650):
    req = requests.post(HOST+'/api/v2/tree-builder/', data=params, verify=False)
    print (req.json())
    try:
        task_id = req.json()['task_id']
    except KeyError:
        return {} , {}
    results = requests.get(HOST + '/api/v2/celery/task/{}'.format(task_id), verify = False)
    clock = 0
    while (not results.json()['complete'] and not results.json().pop('failed', False) and clock <= timeout):
        time.sleep(sleep_time)
        results = requests.get(HOST + '/api/v2/celery/task/{}'.format(task_id), verify = False)
        clock += sleep_time
    return req.json(), results.json()

        
def get_paths (smiles_ls, prioritizer_list, host, params, save_to=None, return_results=True):
    results = {smiles:{} for smiles in smiles_ls}
    params_hist = {smiles:{} for smiles in smiles_ls}
    for smiles in smiles_ls:
        params['smiles'] = smiles
        for prioritizer in prioritizer_list:
            print ("Predicting path for %s using %s" % (smiles, str(prioritizer)))
            params['template_sets'] = prioritizer

            params['template_prioritizers'] = prioritizer
            try:
                request, result = post_and_get(host, params)
            except ConnectionError:
                print ("Connection Error from " + host)
            except Exception as e:
                print ('Failed due to', e)

            params_hist[smiles][str(prioritizer)] = request

            results[smiles][str(prioritizer)] = result
            if save_to != None:
                save_results ({smiles:{str(prioritizer): result}}, 
                          smiles, path = save_to, write_type='a')

        if return_results:
	        return results, params_hist  


def submit_parallel_job (targets, prioritizers, hosts, params, save_prefix, return_results=False):
	n_jobs = len(hosts)
	if return_results:
		parallel_results = Parallel(n_jobs=n_jobs)(delayed(get_paths)(target, prioritizers, host, params, "{}_{}.json".format(save_prefix, i)) for target, host, i in zip(np.array_split(targets, n_jobs), hosts, range(n_jobs)))
		return parallel_results
	else:
		Parallel(n_jobs=n_jobs)(delayed(get_paths)(target, prioritizers, host, params, "{}_{}.json".format(save_prefix, i), return_results=return_results) for target, host, i in zip(np.array_split(targets, n_jobs), hosts, range(n_jobs)))


