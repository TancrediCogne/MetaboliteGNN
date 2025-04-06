# -------------------------------------------- IMPORTS --------------------------------------------
import numpy as np
import pandas as pd
from tqdm import tqdm
import copy
import os

import torch
from torch_geometric.loader import DataLoader
from torch.utils.data import Dataset
import torch.optim as optim
from torch_geometric.nn import GATConv, global_add_pool, GCNConv, GINConv
import torch.nn.functional as F
import torch.nn as nn
from torch.nn import Sequential, Linear, ReLU, BatchNorm1d as BN

from sklearn.metrics import f1_score, recall_score, average_precision_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import matplotlib.pyplot as plt
import seaborn as sns

from datasets import *
from models import *

import warnings
warnings.filterwarnings("ignore")

# -------------------------------------------- TRAINING --------------------------------------------
fingerprint_length = 1024

# TRAIN FUNCTION
def cat_or_start(all, curr, dim=0):
    """Helper function to concatenate tensors or create an empty one

    Args:
        all: all past tensors
        curr: tensor to concatenate to 'all'
        dim: axis along which to concatenatee
    Return: 
        all: concatenated tensors
    """
    if all is None:
        all = curr
    else:
        all = torch.cat((all, curr), dim=dim)
    return all

def train(model, trainloader, validloader, device, criterion, optimizer, threshold, num_epochs, architecture, node_to_keep_idx = None):
    """Training function

    Args:
        model: model used to test 
        trainloader: loader for the training set
        validloader: loader for the validation set
        device: one of ['cuda', 'cpu']
        criterion: loss
        optimizer: Adam optimizer
        threshold: classification threshold
        num_epochs: number of epochs to train the model for
        architecture: one of ['GCN', 'GIN', 'GAT', 'baseline', 'MLP', 'interpretability]
        node_to_keep_idx: None if arch!='interpretability; index of single-output node
    Return: 
        best_epoch: epoch at which the best macro F1-score was obtained
        best_model_wts: weights of the model that yielded the best results
        train_recalls: array containing recall on the training set across epochs
        train_f1_macros: array containing macro F1-score on the training set across epochs
        train_f1_weighteds: array containing weighted F1-score on the training set across epochs
        train_losses: array containing loss on the training set across epochs
        train_aps: array containing area under the precision-recall curve on the training set across epochs
        valid_recalls: array containing recall on the validation set across epochs
        valid_f1_macros: array containing macro F1-score on the validation set across epochs
        valid_f1_weighteds: array containing weighted F1-score on the validation set across epochs
        valid_losses: array containing loss on the validation set across epochs
        valid_aps: array containing area under the precision-recall curve on the validation set across epochs

    """
    model.to(device)

    best_f1_validation = 0.0
    best_epoch = 1
    all_valid_pred = None
    all_valid_true = None
    best_model_wts = copy.deepcopy(model.state_dict())

    train_f1_macros = []
    train_f1_weighteds = []
    train_recalls = []
    train_losses = []
    train_aps = []

    valid_f1_macros = []
    valid_f1_weighteds = []
    valid_recalls = []
    valid_losses = []
    valid_aps = []

    with tqdm(total=num_epochs) as pbar:
        for epoch in range(num_epochs):
            # Training
            model.train()
            train_loss = 0
            y_pred_train = None
            y_true_train = None
            
            for train_data in trainloader:
                optimizer.zero_grad()

                if architecture == 'baseline':
                    outputs = model(train_data[0].to(device))
                    curr_true_y = torch.FloatTensor(np.array(train_data[1].to(device), 'int8'))
                elif architecture == 'MLP':
                    train_x, curr_true_y = train_data[0].to(device), train_data[1].to(device)
                    outputs = model(train_x).squeeze(1)
                elif architecture == 'interpretability':
                    train_data.to(device)
                    outputs = model(train_data.x, train_data.edge_index, train_data.batch, train_data.edge_attr.squeeze(), train_data.graph_embedding).squeeze()
                    curr_true_y = torch.FloatTensor(np.array(train_data.y, 'int8')[:, node_to_keep_idx])
                else:
                    train_data.to(device)
                    outputs = model(train_data.x, train_data.edge_index, train_data.batch, train_data.edge_attr.squeeze(), train_data.graph_embedding)
                    curr_true_y = torch.FloatTensor(np.array(train_data.y, 'int8'))

                loss = criterion(outputs, curr_true_y)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
                curr_outputs = outputs.detach().cpu()

                y_true_train = cat_or_start(y_true_train, curr_true_y)
                y_pred_train = cat_or_start(y_pred_train, curr_outputs)

            # Compute metrics on training set
            y_true_train = y_true_train.cpu()
            y_pred_binary_train = torch.where(torch.sigmoid(y_pred_train) < threshold, torch.tensor(0.0), torch.tensor(1.0))
            train_f1_macro = f1_score(y_true_train, y_pred_binary_train, average="macro")
            train_f1_weighted = f1_score(y_true_train, y_pred_binary_train, average="weighted")
            train_recall = recall_score(y_true_train, y_pred_binary_train, average="macro")
            train_ap = average_precision_score(y_true_train, torch.sigmoid(y_pred_train), average="macro")
            train_f1_macros.append(train_f1_macro)
            train_f1_weighteds.append(train_f1_weighted)
            train_recalls.append(train_recall)
            train_aps.append(train_ap)
            train_losses.append(train_loss / len(trainloader))
    
            pbar.update(1)

            # Validation
            valid_loss = 0.0
            y_true_valid = None
            y_pred_valid = None
            with torch.no_grad():
                model.eval()

                for valid_data in validloader:

                    if architecture == 'baseline':
                        outputs = model(valid_data[0].to(device))
                        curr_true_y = torch.FloatTensor(np.array(valid_data[1].to(device), 'int8'))
                    elif architecture == 'MLP':
                        valid_x, curr_true_y = valid_data[0].to(device), valid_data[1].to(device)
                        outputs = model(valid_x).squeeze(1)
                    elif architecture == 'interpretability':
                        valid_data.to(device)
                        outputs = model(valid_data.x, valid_data.edge_index, valid_data.batch, valid_data.edge_attr.squeeze(), valid_data.graph_embedding).squeeze()
                        curr_true_y = torch.FloatTensor(np.array(valid_data.y, 'int8')[:, node_to_keep_idx])
                    else:
                        valid_data.to(device)
                        outputs = model(valid_data.x, valid_data.edge_index, valid_data.batch, valid_data.edge_attr.squeeze(), valid_data.graph_embedding)
                        curr_true_y = torch.FloatTensor(np.array(valid_data.y, 'int8'))
                    
                    loss = criterion(outputs, curr_true_y)
                    valid_loss += loss.item()
                    curr_outputs = outputs.detach().cpu()
                    
                    y_true_valid = cat_or_start(y_true_valid, curr_true_y)
                    y_pred_valid = cat_or_start(y_pred_valid, outputs)

            # Compute metrics on validation set
            y_true_valid = y_true_valid.cpu()
            y_pred_binary_valid = torch.where(torch.sigmoid(y_pred_valid) < threshold, torch.tensor(0.0), torch.tensor(1.0))

            valid_f1_macro = f1_score(y_true_valid, y_pred_binary_valid, average='macro')
            valid_f1_weighted = f1_score(y_true_valid, y_pred_binary_valid, average='weighted')
            valid_recall = recall_score(y_true_valid, y_pred_binary_valid, average='macro')
            valid_ap = average_precision_score(y_true_valid, torch.sigmoid(y_pred_valid), average="macro")
            valid_f1_macros.append(valid_f1_macro)
            valid_f1_weighteds.append(valid_f1_weighted)
            valid_recalls.append(valid_recall)
            valid_aps.append(valid_ap)
            valid_losses.append(valid_loss / len(validloader))
            
            if all_valid_pred is None:
                all_valid_pred = y_pred_valid
            else:
                all_valid_pred = np.concatenate((all_valid_pred, y_pred_valid))
            if all_valid_true is None:
                all_valid_true = y_true_valid
            else:
                all_valid_true = np.concatenate((all_valid_true, y_true_valid))

            # Update best model weights
            if valid_f1_macro > best_f1_validation:
                    best_f1_validation = valid_f1_macro
                    best_epoch = epoch
                    best_model_wts = copy.deepcopy(model.state_dict())

    return best_epoch, best_model_wts, train_recalls, train_f1_macros, train_f1_weighteds, train_losses, train_aps, valid_recalls, valid_f1_macros, valid_f1_weighteds, valid_losses, valid_aps

# TEST FUNCTION
def test(model, testloader, device, criterion, threshold, architecture, node_to_keep_idx = None):
    """Test function

    Args:
        model: model used to test 
        testloader: loader for the test set
        device: one of ['cuda', 'cpu']
        criterion: loss
        threshold: classification threshold
        architecture: one of ['GCN', 'GIN', 'GAT', 'baseline', 'MLP', 'interpretability]
        node_to_keep_idx: None if arch!='interpretability; index of single-output node
    Return: 
        test_recall: recall on the test set
        test_macro_f1: macro F1-score on the test set
        test_weighted_f1: weighted F1-score on the test set
        test_ap: area under precision-recall curve on the test set
    """
    test_loss = 0.0
    with torch.no_grad():
        model.eval()
        y_true_test = None
        y_pred_test = None

        for test_data in testloader:

            if architecture == 'baseline':
                outputs = model(test_data[0].to(device))
                curr_true_y = torch.FloatTensor(np.array(test_data[1].to(device), 'int8'))
            elif architecture == 'MLP':
                test_x, curr_true_y = test_data[0].to(device), test_data[1].to(device)
                outputs = model(test_x).squeeze(1)
            elif architecture == 'interpretability':
                test_data.to(device)
                outputs = model(test_data.x, test_data.edge_index, test_data.batch, test_data.edge_attr.squeeze(), test_data.graph_embedding).squeeze()
                curr_true_y = torch.FloatTensor(np.array(test_data.y, 'int8')[:, node_to_keep_idx])
            else:
                test_data.to(device)
                outputs = model(test_data.x, test_data.edge_index, test_data.batch, test_data.edge_attr.squeeze(), test_data.graph_embedding)
                curr_true_y = torch.FloatTensor(np.array(test_data.y, 'int8'))

            loss = criterion(outputs, curr_true_y)
            test_loss += loss.item()
            curr_outputs = outputs.detach().cpu()

            y_true_test = cat_or_start(y_true_test, curr_true_y)
            y_pred_test = cat_or_start(y_pred_test, curr_outputs)
    
        # Compute metrics on test set
        y_true_test = y_true_test.cpu()
        y_pred_binary_test = torch.where(torch.sigmoid(y_pred_test) < threshold, torch.tensor(0.0), torch.tensor(1.0))
        test_macro_f1 = f1_score(y_true_test, y_pred_binary_test, average='macro')
        test_weighted_f1 = f1_score(y_true_test, y_pred_binary_test, average='weighted')
        test_recall = recall_score(y_true_test, y_pred_binary_test, average='macro')
        test_ap = average_precision_score(y_true_test, torch.sigmoid(y_pred_test), average='macro')

        return test_recall, test_macro_f1, test_weighted_f1, test_ap

# CROSS VALIDATION
def cross_validation(dataset_train, dataset_test, device, threshold, num_epochs, lr, weight_decay, pos_weight, num_folds, architecture, level1_node_num, chemberta, node_to_keep_idx= None):
    """Cross-validation function

    Args:
        dataset_train: training dataset that will be split intro train/valid 
        dataset_test: test set
        device: one of ['cuda', 'cpu']
        threshold: classification threshold
        num_epochs: number of epochs to run the training for
        lr: learning rate
        weight_decay: penalty parameter
        pos_weight: weight of positive examples
        num_folds: number of cross-validation folds
        architecture: one of ['GCN', 'GIN', 'GAT', 'baseline', 'MLP', 'interpretability]
        level1_node_num: one of ['disposition', 'role', 'process', 'physiological_effect]
        chemberta: if using ChemBERTa embeddings (only compatible with GNN models)
        node_to_keep_idx: None if arch!='interpretability; index of single-output node
    Return: 
        bestest_model_wts: best model weights after training
        all_train_recalls: array containing recall of all folds on training set
        all_train_macro_f1s: array containing macro F1-score of all folds on training set
        all_train_weighted_f1s: array containing weighted F1-score of all folds on training set
        all_train_losses: array containing loss of all folds on training set
        all_train_aps: array containing area under precision-recall curve of all folds on training set
        all_valid_recalls: array containing recall of all folds on validation set
        all_valid_macro_f1s: array containing macro F1-score of all folds on validation set
        all_valid_weighted_f1s: array containing weighted F1-score of all folds on validation set
        all_valid_losses: array containing loss of all folds on validation set
        all_valid_aps: array containing area under precision-recall curve of all folds on validation set
        all_test_recalls: array containing recall of all folds on test set
        all_test_macro_f1s: array containing macro F1-score of all folds on test set
        all_test_weighted_f1s: array containing weighted F1-score of all folds on test set
        all_test_aps: array containing area under precision-recall curve of all folds on test set
    """
    skf = KFold(n_splits=num_folds, shuffle=True, random_state=45) # Compute folds
    batch_size = 64

    testloader = DataLoader(dataset_test, batch_size=batch_size) # Load test set

    all_train_macro_f1s = np.zeros((num_folds, num_epochs))
    all_train_weighted_f1s = np.zeros((num_folds, num_epochs))
    all_train_recalls = np.zeros((num_folds, num_epochs))
    all_train_losses = np.zeros((num_folds, num_epochs))
    all_train_aps = np.zeros((num_folds, num_epochs))

    all_valid_macro_f1s = np.zeros((num_folds, num_epochs))
    all_valid_weighted_f1s = np.zeros((num_folds, num_epochs))
    all_valid_recalls = np.zeros((num_folds, num_epochs))
    all_valid_losses = np.zeros((num_folds, num_epochs))
    all_valid_aps = np.zeros((num_folds, num_epochs))

    all_test_macro_f1s = []
    all_test_weighted_f1s = []
    all_test_recalls = []
    all_test_aps = []

    max_f1 = 0
    bestest_model_wts = None

    curr_fold = 0
    # Initialize the model for this run (inside the cross-validation loop)
    for train_index, val_index in skf.split(dataset_train):
        print(f"***************************************************** Fold {curr_fold+1}/{num_folds} *****************************************************")
        if architecture == 'GCN':
            model = GCN(num_node_features=63, hidden_dim1=32, hidden_dim2=32, output_dim=level1_node_num, chemberta=chemberta) 
        elif architecture == 'GIN':
            model = GIN(num_features = 63, num_layers=2, hidden_dim=32, output_dim=level1_node_num, chemberta=chemberta)
        elif architecture == 'GAT':
            model = GAT(num_features = 63, heads=8, hidden_channels=32, output_dim=level1_node_num, chemberta=chemberta)
        elif architecture == 'MLP':
            model = MLP(input_dim=graph_embedding_length, hidden_dim1=32, hidden_dim2=32, output_dim=level1_node_num)
        elif architecture == 'baseline':
            model = MLP(input_dim=fingerprint_length, hidden_dim1=32, hidden_dim2=32, output_dim=level1_node_num)
        elif architecture == 'interpretability':
            model = GAT(num_features = 63, heads=8, hidden_channels=32, output_dim=1, chemberta=chemberta)
            
        criterion = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        criterion.to(device)
        optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

        # Split training set into training and validation sets for each fold
        fold_train_dataset = [dataset_train[i] for i in train_index]

        fold_valid_dataset = [dataset_train[i] for i in val_index]
        fold_trainloader = DataLoader(fold_train_dataset, batch_size=batch_size, shuffle=True)

        fold_validloader = DataLoader(fold_valid_dataset, batch_size=batch_size)

        # Train and evaluate the model for this fold
        if architecture == 'MLP':
            results = train(model, ChemBERTA_Dataset(fold_trainloader), ChemBERTA_Dataset(fold_validloader), device, criterion, optimizer, threshold, num_epochs, architecture)
        else:
            results = train(model, fold_trainloader, fold_validloader, device, criterion, optimizer, threshold, num_epochs, architecture, node_to_keep_idx=node_to_keep_idx)
        
        _, best_model_wts, train_recalls, train_macro_f1s, train_weighted_f1s, train_losses, train_aps, valid_recalls, valid_macro_f1s, valid_weighted_f1s, valid_losses, valid_aps = results

        # Append current fold's metrics
        all_train_macro_f1s[curr_fold] = train_macro_f1s
        all_train_weighted_f1s[curr_fold] = train_weighted_f1s
        all_train_recalls[curr_fold] = train_recalls
        all_train_losses[curr_fold] = train_losses
        all_train_aps[curr_fold] = train_aps

        all_valid_macro_f1s[curr_fold] = valid_macro_f1s
        all_valid_weighted_f1s[curr_fold] = valid_weighted_f1s
        all_valid_recalls[curr_fold] = valid_recalls
        all_valid_losses[curr_fold] = valid_losses
        all_valid_aps[curr_fold] = valid_aps

        # Load best model weights and test on the current fold
        model.load_state_dict(best_model_wts)
        if architecture == 'MLP':
            test_recall, test_macro_f1, test_weighted_f1, test_ap = test(model, ChemBERTA_Dataset(testloader), device, criterion, threshold, architecture)
        else:
            test_recall, test_macro_f1, test_weighted_f1, test_ap = test(model, testloader, device, criterion, threshold, architecture, node_to_keep_idx=node_to_keep_idx)
        
        # Append test metrics
        all_test_recalls.append(test_recall)
        all_test_macro_f1s.append(test_macro_f1)
        all_test_weighted_f1s.append(test_weighted_f1)
        all_test_aps.append(test_ap)

        # Update best model weights
        if test_macro_f1 > max_f1:
            max_f1 = test_macro_f1
            bestest_model_wts = best_model_wts
        curr_fold+=1

    return bestest_model_wts, all_train_recalls, all_train_macro_f1s, all_train_weighted_f1s, all_train_losses, all_train_aps, all_valid_recalls, all_valid_macro_f1s, all_valid_weighted_f1s, all_valid_losses, all_valid_aps, all_test_recalls, all_test_macro_f1s, all_test_weighted_f1s, all_test_aps

# PLOT METRIC
def plot_metric(train_values, valid_values, metric_name, saving_folder, log=False):
    """Plot a given metric across epochs

    Args:
        train_values: values of the given metric on the training set
        valid_values: values of the given metric on the validation set
        metric_name: name of the metric ti plot
        saving_foler: folder where to store the plot
        log: boolean to indicate if y-axis must be log-scaled
    """
    mean_train_values = np.mean(train_values, axis=0)
    mean_val_values = np.mean(valid_values, axis=0)
    palette = sns.color_palette("mako_r", 2)
    epochs = np.arange(1, len(mean_train_values) + 1)

    plt.figure(figsize=(10, 6))

    sns.lineplot(x=np.concatenate([epochs, epochs]),
                 y=np.concatenate([mean_train_values, mean_val_values]),
                 hue=['Training'] * len(mean_train_values) + ['Validation'] * len(mean_val_values),
                 palette = palette,
                 errorbar='sd')

    plt.xlabel('Epochs')
    plt.xticks(range(0, len(mean_train_values) + 1, 5))
    plt.grid(True, linestyle='--', alpha=0.7)
    if log:
        plt.yscale('log')
        plt.ylabel(f'{metric_name} (log-scale)')
    else:
        plt.ylabel(metric_name)

    plt.legend()
    plt.savefig(saving_folder + f"/evolution_{metric_name.lower()}")
    plt.close()

# MODEL TRAINING
def model_training(dataset_date, saving_folder, level1_node, level1_node_num, architecture, chemberta, threshold, num_epochs, lr, weight_decay, node_to_keep_idx=None):
    """Train model, save best parameters and plots, and returns metrics

    Args:
        dataset_date: date of data processing
        saving_folder: general folder where all the models will be stored
        level1_node: one of ['disposition', 'role', 'process', 'physiological_effect]
        level1_node_num: 31 (disposition), 11 (role), 14 (process), 16 (physiological effect)
        architecture: one of ['GCN', 'GIN', 'GAT', 'baseline', 'MLP', 'interpretability]
        chemberta: if using ChemBERTa embeddings (only compatible with GNN models)
        threshold: classification threshold
        num_epochs: number of epochs to run the training for
        lr: learning rate
        weight_decay: penalty parameter
        node_to_keep_idx: None if arch!='interpretability; index of single-output node 
    Return: 
        average_recall: average recall across 5 validation folds
        std_recall: standard deviation of recall across 5 validation folds
        average_f1_macro: average macro F1-score across 5 validation folds
        f1_macro_std: standard deviation of macro F1-score across 5 validation folds
        average_f1_weighted: average weighted F1-score across 5 validation folds
        f1_weighted_std: standard deviation of weighted F1-score across 5 validation folds
        average_ap: average area under precision-recall curve across 5 validation folds
        std_ap: standard deviation of area under precision-recall curve across 5 validation folds
    """
    os.makedirs(saving_folder, exist_ok=True) # Create folder where to store model results 
    # Initialize the model, loss function, and optimizer
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load datasets
    if architecture == 'baseline':
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=fingerprint_length)
        dataset_train, dataset_test = data_generation_fingerprints(level1_node, dataset_date, morgan_gen)
    else:   
        dataset_train = torch.load('data/processed_data/processed_data_' + dataset_date + '_' + level1_node +'/dataset_train.pt')
        dataset_test = torch.load('data/processed_data/processed_data_' + dataset_date + '_' + level1_node +'/dataset_test.pt')
    pos_weight = torch.load('data/processed_data/processed_data_' + dataset_date + '_' + level1_node +'/pos_weight.pt')
    
    results = cross_validation(dataset_train, dataset_test, device, threshold, num_epochs,  lr, weight_decay, pos_weight, num_folds=5, architecture=architecture, level1_node_num=level1_node_num, chemberta=chemberta, node_to_keep_idx=node_to_keep_idx)
    best_model_wts, all_train_recalls, all_train_f1_macros, all_train_f1_weighteds, all_train_losses, all_train_aps, all_valid_recalls, all_valid_f1_macros, all_valid_f1_weighteds, all_valid_losses, all_valid_aps, all_test_recalls, all_test_f1_macros, all_test_f1_weighteds, all_test_aps  = results

    # Save best model weights 
    torch.save(best_model_wts, saving_folder +"/best_model_wts.pt")

    # Plot metrics
    plot_metric(all_train_recalls, all_valid_recalls, "Recall", saving_folder, log=False)
    plot_metric(all_train_f1_macros, all_valid_f1_macros, "Macro F1-Score", saving_folder, log=False)
    plot_metric(all_train_f1_weighteds, all_valid_f1_weighteds, "Weighted F1-Score", saving_folder, log=False)
    plot_metric(all_train_aps, all_valid_aps, "Average Precision", saving_folder, log=False)
    plot_metric(all_train_losses, all_valid_losses, "Loss", saving_folder, log=True)

    # Print final cross-validation results
    average_recall = np.mean(all_test_recalls)
    recall_std = np.std(np.array(all_test_recalls))
    average_f1_macro = np.mean(all_test_f1_macros)
    f1_macro_std = np.std(np.array(all_test_f1_macros))
    average_f1_weighted = np.mean(all_test_f1_weighteds)
    f1_weighted_std = np.std(np.array(all_test_f1_weighteds))
    average_ap = np.mean(all_test_aps)
    ap_std = np.std(np.array(all_test_aps))

    return average_recall, recall_std, average_f1_macro, f1_macro_std, average_f1_weighted, f1_weighted_std, average_ap, ap_std
