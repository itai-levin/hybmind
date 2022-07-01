This repository contains a collection of scripts and jupyter notebooks used to perform data processing and analysis for the manuscript "Merging enzymatic and synthetic chemistry with computational synthesis planning." This repository does not include code for to perform synthesis planning which can be found at 

## Installation
There is no installation required. The scripts can be run independently after the repository is cloned and the necessary python packages are installed.  

## Overview
The jupyter notebooks can be run by clicking through them. Example commands are provided for the bash and python scripts. All code was run locally unless otherwise noted.

1. Process BKMS database:
* `process_BKMS_database/Notebooks/process_bkms_db.ipynb`: extracts unmapped reaction SMILES and performs initial removal of cofactors. Generates a .csv file with the input BKMS data file and an additional column with unmapped reaction SMILES strings. Estimated runtime: 5 minutes. 
* `process_BKMS_database/scripts/parallel_map.sh`: performs atom-atom mapping. Generates .csv files of a portion of the processed BKMS data file with an additional column containing atom-atom mapped reaction SMILES strings. The number of files generated is the same as the number of tasks set. To run: `bash parallel_map.sh -i ../../data/processed_bkms_database/bkms_w_SMILES.txt -o ../../data/processed_bkms_database/bkms -t 4 -c 21`. Run remotely. Estimated runtime: 3 hours.
* `process_BKMS_database/scripts/combine_mapped_files.sh`: aggregates the atom-atom mapped reactions. To run: `bash combine_mapped_files.sh -d ../../data/processed_bkms_database/ -p bkms- -o mapped -x ../../data/processed_bkms_database/`. Estimated runtime: 5 seconds.
* `process_BKMS_database/Notebooks/process_mapped_bkms.ipynb`: performs an additional round of co-factor removal to yield single-product reactions. Estimated runtime: 1 hour.

2. Training template prioritizer:
More detailed documention of functions here[https://gitlab.com/mefortunato/template-relevance]
* `train-model/hyperparameter_optimization.sh`: runs grid search to train models. To run: `bash hyperparameter_optimization.sh -p bkms- -m bkms -t bkms`. Run remotely. Estimated runtime: 2 hours. 
* `process.py` and `train.py` used to train the final model on all of the processed BKMS data. Estimated runtime: 1 hour.

From here, trained models were copied to a GCP VM instance running a deployment of ASKCOS [(separate repo)](https://github.com/itai-levin/chemoenzymatic-askcos)

3. Analysis of templates:
*  `analyze_templates/scripts/count_templates.ipynb`: counts the number of precedents for each template and the effect of setting precedent thresholds for the templates. Data presented in Figure 2. Estimated runtime: 6 hours.
* `analyze_templates/scripts/embarassingly_parallel_template_application.sh`: applies all templates extracted from Reaxys to the reactions from BKMS and returns which reactions can be recovered. To run: `bash embarassingly_parallel_template_application.sh`. Run remotely. Estimated runtime: 20 hours.
* `analyze_templates/Notebooks/unique enzyme chemistry.ipynb`. Analysis of the outputs from the template application script. Data is partially presented in Figure 3. Estimated runtime: 10 minutes.
* `analyze_templates/Notebooks/compare_model_outputs.ipynb`: compares the magnitude of the scores output from the model trained on Reaxys and the model trained on BKMS by inputting molecule sets from MOSES and ZINC. This data is presented in Figure 4. Estimated runtime: 2 hours.

4. Calling pathway searches from the API
* call_ASKCOS_api/submit_tree_search.py: submits requests to deployed ASKCOS servers to generate and store tree search graphs. To run: `python submit_tree_search.py`. Run remotely. Estimated runtime: 50 hours.

5. Processing the results of the ASKCOS tree builder search
* process_ASKCOS_results/Notebooks/parse_returned_synthesis_graphs_boutiques: extracts the information from the synthesis routes (available in compressed form at `data/boutique_1000_retrosynthesis_graphs.json.gz`) presented in Figure 6. Run remotely. Estimated runtime: 30 minutes.

The final templates, reactions, and trained model used in our ASKCOS searches are available in a [separate repository](https://github.com/itai-levin/bkms-data)

### System information

For code that was run locally:
* MacOS 10.14.6 (Mojave) with Python 3.7.9
* 8 GB RAM, 1.7GHz Intel Core i7, 2 cores

For code that was run remotely:
* Ubuntu 20.04.1 LTS with Python 3.9.7
* 126 GB RAM, 64 cores

#### Python libraries used
  - rdkit: 2020.03.2
  - pandas: 1.0.5
  - matplotlib: 3.1.3
  - numpy: 1.19.5
  - tensorflow: 2.1.0
  - requests: 2.27.1


