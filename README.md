## Data pre-processing
- File needed: `hmdb_metabolites.xml` (can be downloaded at: https://hmdb.ca/downloads)
1. Split the xml file into multiple files by running: `python3 xml_breaker.py hmdb_metabolites.xml tag number` 
    - `hmdb_metabolites.xml`
    - `tag`: the tag name where the splits will be made (in this case: metabolite)
    - `number`: the number of such tags per output file (in this case: around 2500 should work, but it might depend on the OS)
2. Extract the tree and truth table of each smaller xml file by running: `python3 extract_trees_and_truth_tables.py hmdb_metabolitesX.xml truth_tables trees`
    - `hmdb_metabolitesX.xml`: the xml file from which the tree and truth table will be extracted (X represents the number of the file)
    - `truth_tables`: the path to the folder where all the truth tables will be saved
    - `trees`: the path to the folder where all the trees will be saved
3. Merge all the trees and truth tables and extract all infos by running: `python3 extract_all_infos.py hmdb_metabolites.xml trees truth_tables`:
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
- Process the data by running: `python3 data_processing.py structures.sdf filtered_disposition_truth_table.csv data_folder ids.npy filtered_nodes_mad.npy smiles_quantified.csv`
    - `structures.sdf`: file containing all infos to build graphs (can be downloaded at: https://hmdb.ca/downloads)
    - `filtered_disposition_truth_table.csv`: file containing the truth table for given level 1 and "detected and quantified" ids
    - `data_folder`: folder where to store processed data
    - `quantified_ids.npy`: "detected and quantified" ids
    - `filtered_nodes_mad.npy`: nodes selected with MAD filtering
    - `smiles_quantified.csv`: smiles of each metabolite in detected and quantified category, used to retrieve chemBERTa embeddings
    
- This procedure will result in six files stored in the folder `processed_data_date_level_1_node`:
    - `adjacency_features.npy`: array of shape [G, 2, N] containing information on edges (from-to) (G: # of graphs)
    - `dataset.pt`: dataset saved as a pytorch geometric object
    - `edge_features.npy`: array of shape [G, N] containing info on type of bond
    - `filtered_outputs.csv`: contains the y_true values containing the exact same information stored in dataset.y
    - `node_features.npy`: array of shape [G, M, P] containing 2D coordinnates + one-hot encoded atomic number for each atom
    - `valid_ids.npy`: contains all the ids of the metabolites used in `dataset.pt`

## Model training
- Files needed:
    - `dataset.pt`
    - `filtered_outputs.csv`
    
- This procedure will result in five files stored in the folder `model_date_level_1_node`:
    - `best_model_wts.pt`: the best model weights to use for interpretability
    - `evolution_accuracy.png`: plot across epochs to show evolution of accuracy
    - `evolution_f1-score.png`: plot across epochs to show evolution of f1-score
    - `evolution_loss.png`: plot across epochs to show evolution of loss
    - `test_idx.npy`: list of idx of each metabolite in the test set
 
## Interpretability
- 
