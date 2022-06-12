## Overview
1. Process BKMS database:
* `process_BKMS_database/Notebooks/process_bkms_db.ipynb`: extracts unmapped reaction SMILES and performs initial removal of cofactors
* `process_BKMS_database/scripts/parallel_map.sh`: performs atom-atom mapping 
* `process_BKMS_database/scripts/combine_mapped_files.sh`: aggregates the atom-atom mapped reactions
* `process_BKMS_database/Notebooks/process_mapped_bkms.ipynb`: performs an additional round of co-factor removal to yield single-product reactions

2. Training template prioritizer:
More detailed documention of functions here[https://gitlab.com/mefortunato/template-relevance]
* `train-model/hyperparameter_optimization.sh`: runs grid search to train models
* `process.py` and `train.py` used to train the final model on all of the processed BKMS data

From here, trained models were copied to a GCP VM instance running a deployment of ASKCOS [(separate repo)](https://github.com/itai-levin/chemoenzymatic-askcos)

3. Analysis of templates:
*  `analyze_templates/scripts/count_templates.ipynb`: counts the number of precedents for each template and the effect of setting precedent thresholds for the templates. Data presented in Figure 2. 
* `analyze_templates/scripts/embarassingly_parallel_template_application.sh`: applies all templates extracted from Reaxys to the reactions from BKMS and returns which reactions can be recovered. Analysis of the outputs is done in `analyze_templates/Notebooks/unique enzyme chemistry.ipynb`. Data is partially presented in Figure 3.
* `analyze_templates/Notebooks/compare_model_outputs.ipynb`: compares the magnitude of the scores output from the model trained on Reaxys and the model trained on BKMS by inputting molecule sets from MOSES and ZINC. This data is presented in Figure 4.

4. Calling pathway searches from the API
* call_ASKCOS_api/submit_tree_search.py: submits requests to deployed ASKCOS servers to generate and store tree search graphs.

5. Processing the results of the ASKCOS tree builder search
* process_ASKCOS_results/Notebooks/parse_returned_synthesis_graphs_boutiques: extracts the information from the synthesis routes (available in compressed form at `data/boutique_1000_retrosynthesis_graphs.json.gz` ) presented in Figure 6.

The final templates, reactions, and trained model used in our ASKCOS searches are available in a [separate repository](https://github.com/itai-levin/bkms-data)
