# IMPORTS
import numpy as np
import pandas as pd
from tqdm import tqdm
import copy
import torch
from torch.utils.data import random_split, Dataset
from sklearn.preprocessing import StandardScaler

# CLASS DEF
class GraphDataset(Dataset):
    def __init__(self, data_list):
        self.dataset = data_list

    def __len__(self):
        return len(self.dataset)

    def __getitem__(self, idx):
        return self.dataset[idx]
    def get_dataset(self):
        return copy.deepcopy(self.dataset)

# STRATIFY
def stratify(data, idx_dict, ratio):
    """Stratify data based on grouped metabolites

    Args:
        data: data to split
        idx_dict: dictionary given a list of metabolite index for each group
        ratio: ratio between train and test
    Return: 
        train_data: training data
        test_data: test data
    """
    test_data = []
    train_data = []
    rest_data = []

    # Iterate over all groups and split each into train/test
    for key, val in idx_dict.items():
        if len(val) >= 10:
            curr_output = [data[v] for v in val.values]
            curr_train, curr_test = random_split(curr_output, ratio)
            train_data += curr_train
            test_data += curr_test
        # Groups that have less than 10 metabolites can not be split properly and are all added to a bigger group before a final split
        else:
           rest_data += [data[v] for v in val.values]

    # Split the group containing all 'under 10 metabolites' groups
    curr_train, curr_test = random_split(rest_data, ratio)
    train_data += curr_train
    test_data += curr_test
    train_to_remove = []
    test_to_remove = []
    
    # Remove all metabolites that have only one bond
    for t in train_data:
        if t.edge_attr.shape[0] == 1:
            train_to_remove.append(t)
    for t in test_data:
        if t.edge_attr.shape[0] == 1:
            test_to_remove.append(t)

    train_data = [t for t in train_data if t not in train_to_remove]
    test_data = [t for t in test_data if t not in test_to_remove]
    return train_data, test_data

# SPLIT SETS
def output_based_split():
    """Split dataset into training and test seets

        Return: 
            dataset_train: training set
            dataset_test: test set
            pos_weight: weight of positive examples
        """
    # Load data and scale it
    dataset = torch.load('data/processed_data/processed_data_' + dataset_date + '_' + level1_node +'/dataset.pt')
    for d in dataset:
        x = d.x[:, :2]
        scaler = StandardScaler()
        scaler.fit(x)
        scaled_x = scaler.transform(x)
        d.x[:,:2] = torch.tensor(scaled_x)

    # Load and group metabolites based on outputs
    y_true = pd.read_csv('data/processed_data/processed_data_' + dataset_date + '_' + level1_node +'/filtered_outputs.csv').fillna(False)
    df = y_true.drop(['accession'], axis=1)
    df['RowString'] = df.astype(str).agg(','.join, axis=1)
    index_dict = df.groupby('RowString').groups

    # Stratify data based on groups 
    dataset_train, dataset_test = stratify(dataset, index_dict, [0.9,0.1]) 

    # Calculatee weight of positive examples
    train_ys = [d.y for d in dataset_train]
    train_ys_df = pd.DataFrame(train_ys)
    num_pos = train_ys_df.value_counts().values[0]
    num_neg = len(train_ys_df) - num_pos
    pos_weight = num_neg / num_pos
    pos_weight = torch.tensor([pos_weight], dtype=torch.float)
    return dataset_train, dataset_test, pos_weight

dataset_date = "25_03_31" # TOCHANGE
level1_node = "disposition" # one of ['disposition', 'role', 'process', 'physiological_effect']
level1_node_num = 31 # 31 (disposition), 11 (role), 14 (process), 16 (physiological effect)

# Load dataset before splitting
dataset = torch.load('data/processed_data/processed_data_' + dataset_date + '_' + level1_node +'/dataset.pt')

dataset_train, dataset_test, pos_weight = output_based_split() # Split dataset into train/test

print("Saving training dataset of size: ", len(dataset_train))
torch.save(dataset_train, 'data/processed_data/processed_data_' + dataset_date + '_' + level1_node + "/dataset_train.pt")
print("Saving test dataset of size: ", len(dataset_test))
torch.save(dataset_test, 'data/processed_data/processed_data_' + dataset_date + '_' + level1_node + "/dataset_test.pt")
torch.save(pos_weight, 'data/processed_data/processed_data_' + dataset_date + '_' + level1_node + "/pos_weight.pt")
