# MetaboliteGNN 
This repository contains all the code asssociated with the project (TODO link). It allows to train multi-label classifier for metabolite function prediction based on structural information. All the processing and training steps are explained below.

## Data pre-processing
- File needed: `hmdb_metabolites.xml` (can be downloaded at: https://hmdb.ca/downloads)
1. Split the xml file into multiple files by running:
   ```
   python3 xml_breaker.py hmdb_metabolites.xml tag number
   ```
    - `hmdb_metabolites.xml`
    - `tag`: the tag name where the splits will be made (in this case: metabolite)
    - `number`: the number of such tags per output file (in this case: around 2500 should work, but it might depend on the OS)
3. Extract the tree and truth table of each smaller xml file by running:
   ```
   python3 extract_trees_and_truth_tables.py hmdb_metabolitesX.xml truth_tables trees
   ```
    - `hmdb_metabolitesX.xml`: the xml file from which the tree and truth table will be extracted (X represents the number of the file)
    - `truth_tables`: the path to the folder where all the truth tables will be saved
    - `trees`: the path to the folder where all the trees will be saved
4. Merge all the trees and truth tables and extract all infos by running: `python3 extract_all_infos.py hmdb_metabolites.xml trees truth_tables`:
    - `hmdb_metabolites.xml`: the xml file containing the whole database
    - `truth_tables`: the path to the folder where all the truth tables will be saved
    - `trees`: the path to the folder where all the trees will be saved

- This procedure will result in three files stored in the folder `preprocessed_data`:
    - `ontology_tree.csv`
    - `ontology_truth_table.csv`
    - `metabolites_infos.csv`

## Data processing
- Files needed:
    - `ontology_tree.csv`
    - `ontology_truth_table.csv`
    - `metabolites_infos.csv`
    - `structures.sdf` (can be downloaded at: https://hmdb.ca/downloads)
- Process the data by running:
  ```
  python3 data_processing.py structures.sdf filtered_disposition_truth_table.csv data_folder ids.npy filtered_nodes_mad.npy smiles_quantified.csv
  ```
    - `structures.sdf`: file containing all infos to build graphs
    - `filtered_disposition_truth_table.csv`: file containing the truth table for given level 1 and "detected and quantified" ids (can be found in `Data` folder)
    - `data_folder`: folder where to store processed data
    - `quantified_ids.npy`: "detected and quantified" ids (can be found in `Data` folder)
    - `filtered_nodes_mad.npy`: nodes selected with MAD filtering (can be found in `Data` folder)
    - `smiles_quantified.csv`: smiles of each metabolite in detected and quantified category, used to retrieve chemBERTa embeddings (can be found in `Data` folder)  
- This procedure will result in six files stored in the folder `processed_data_date_level_1_node`:
    - `adjacency_features.npy`: array of shape [G, 2, N] containing information on edges (from-to) (G: # of graphs)
    - `dataset.pt`: dataset saved as a pytorch geometric object
    - `edge_features.npy`: array of shape [G, N] containing info on type of bond
    - `filtered_outputs.csv`: contains the y_true values containing the exact same information stored in dataset.y
    - `node_features.npy`: array of shape [G, M, P] containing 2D coordinnates + one-hot encoded atomic number for each atom
    - `valid_ids.npy`: contains all the ids of the metabolites used in `dataset.pt`

## Train/test split
- Files needed:
    - `dataset.pt`: dataset saved as a pytorch geometric object
 - Parameters to choose:
    - `dataset_date`: date at which the data was processed (this assumes that the processed data is in a folder named `processed_data_dataset_date`)
    - `level1_node`: one of 'disposition', 'role', 'process', 'physiological_effect' (based on which data is being processed)
    - `level1_node_num`: depends on which `level_1_node` was chosen (31 (disposition), 11 (role), 14 (process), 16 (physiological effect))
- Split the data by running:
  ```
  python3 train_test_split.py
  ```
- This procedure will result in three files:
    - `dataset_train.pt`: training dataset saved as a pytorch geometric object
    - `dataset_test.pt`: test dataset saved as a pytorch geometric object
    - `pos_weight.pt`: weight of positive examples asved as a pytorch geometric object
      
## Model training 
- The file `model.py` is meant to train multiple models with different architectures for comparison. The file `model.ipynb` can be used to train one single model (described below)
- Parameters to choose:
    - `dataset_date`: date at which the data was processed (this assumes that the processed data is in a folder named `processed_data_dataset_date`)
    - `current_date`: date at which the model is trained
    - `level1_node`: one of 'disposition', 'role', 'process', 'physiological_effect' (based on which data is being processed)
    - `level1_node_num`: depends on which `level_1_node` was chosen (31 (disposition), 11 (role), 14 (process), 16 (physiological effect))
    - `chemberta`: whether to use ChemBERTa embeddings or not
    - `node_to_keep`: name of the single-node to predict (only useful if training interpretability model)
    - `node_to_keep_idx`: index of the single-node to predict in the output (only useful if training interpretability model)
    - `arch`: one of 'GCN', 'GIN', 'GAT', 'baseline', 'MLP', 'interpretability' 
- This procedure will result in six files stored in the folder `model_date_level_1_node_arch_chemberta`:
    - `best_model_wts.pt`: the best model weights to use for interpretability
    - `evolution_average precision.png`: plot across epochs to show evolution of area under precision-recall curve
    - `evolution_loss.png`: plot across epochs to show evolution of loss
    - `evolution_macro_f1-score.png`: plot across epochs to show evolution of macro f1-score
    - `evolution_recall.png`: plot across epochs to show evolution of recall
    - `evolution_weighted_f1-score.png`: plot across epochs to show evolution of weighted f1-score
