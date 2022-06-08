## Overview
1. Process BKMS database:
* process_BKMS_database/Notebooks/process_bkms_db.ipynb: extracts unmapped reaction SMILES and performs initial removal of cofactors
* process_BKMS_database/scripts/parallel_map.sh: performs atom-atom mapping 
* process_BKMS_database/scripts/combine_mapped_files.sh: aggregates the atom-atom mapped reactions
* process_BKMS_database/Notebooks/process_mapped_bkms.ipynb: performs an additional round of co-factor removal to yield single-product reactions

2. Training template prioritizer:
More detailed documention of functions here[https://gitlab.com/mefortunato/template-relevance]
* train-model/hyperparameter_optimization.sh: runs grid search to train models

From here, trained models were copied to a GCP VM instance running a deployment of ASKCOS (separate repo)

3. Analysis of templates:
* analyze_templates/scripts/embarassingly_parallel_template_application.sh: applies all templates extracted from Reaxys to the reactions from BKMS and returns which reactions can be recovered.
* 

4. Calling pathway searches from the API
* call_ASKCOS_api/submit_tree_search.py: submits requests to deployed ASKCOS servers to generate and store tree search graphs.

5. Processing the results of the ASKCOS tree builder search
* process_ASKCOS_results/Notebooks/parse_returned_synthesis_graphs_boutiques: extracts the information from the synthesis routes presented in Figure 5.