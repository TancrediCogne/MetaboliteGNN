# -------------------------------------------- IMPORTS --------------------------------------------
import pandas as pd
import os

from rdkit.Chem import rdFingerprintGenerator

from datasets import *
from models import *
from model_helpers import *

import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------- MAIN ----------------------------------------------
# Model parameters
dataset_date = "25_03_31" 
current_date = "25_04_05"
chemberta = False # if using ChemBERTa embeddings (only compatible with GNN models)
level1_nodes = ['disposition', 'role', 'process', 'physiological_effect'] 
level1_node_nums = [31, 11, 14, 16]
architectures = ['GCN', 'GIN', 'GAT'] if chemberta else ['baseline', 'MLP', 'GCN', 'GIN', 'GAT']
os.makedirs("models_" + dataset_date, exist_ok=True)

# Hyper-parameters
num_epochs = 1
threshold = 0.5
lr = 0.005
weight_decay=0.001

# Dataframe to store all the results on test set
all_results = pd.DataFrame(columns=['Name', 'Avg Recall', 'Std Recall', 
                                            'Avg Macro F1-score', 'Std Macro F1-score', 
                                            'Avg weighted F1-score', 'Std weighted F1-score', 
                                            'Avg AP', 'Std AP'])

# Model training/testing 
for i, level1_node in enumerate(level1_nodes):
    print(f'--------------------------------------------------------------------------------------------------------------------')
    print(f'STARTING TO TRAIN "{str.upper(level1_node)}" MODELS')
    for arch in architectures:
        print(f'-----------------------------------------------Architecture set to {arch}-----------------------------------------------')
        if chemberta:
            saving_folder = "models_" + current_date + "/model_" + level1_node + "_" + arch + '_chemberta' 
        else:
            saving_folder = "models_" + current_date + "/model_" + level1_node + "_" + arch
        average_recall, recall_std, average_f1_macro, f1_macro_std, average_f1_weighted, f1_weighted_std, average_ap, ap_std = model_training(dataset_date, saving_folder, level1_node, level1_node_nums[i], arch, chemberta, threshold, num_epochs, lr, weight_decay)
        all_results.loc[len(all_results.index)] = [level1_node + '_' + arch, average_recall, recall_std, average_f1_macro, f1_macro_std, average_f1_weighted, f1_weighted_std, average_ap, ap_std]
        print(f'\n RESULTS: Average Macro F1-score: {average_f1_macro:.3f} ± {f1_macro_std:.3f}')
        print(f'          Average Weighted F1-score: {average_f1_weighted:.3f} ± {f1_weighted_std:.3f}')
        print(f'          Average Recall: {average_recall:.3f} ± {recall_std:.3f}')
        print(f'          Average AP: {average_ap:.3f} ± {ap_std:.3f}\n')

# Storing results in a .csv file
if chemberta:
    all_results.to_csv('models_' + current_date + '/chemberta.csv')
else:
    all_results.to_csv('models_' + current_date + '/no_chemberta.csv')