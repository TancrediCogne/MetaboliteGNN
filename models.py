# -------------------------------------------- IMPORTS --------------------------------------------
import numpy as np

import torch
from torch_geometric.nn import GATConv, global_add_pool, GCNConv, GINConv
import torch.nn.functional as F
import torch.nn as nn
from torch.nn import Sequential, Linear, ReLU, BatchNorm1d as BN

import warnings
warnings.filterwarnings("ignore")

# --------------------------------------------- MODELS ---------------------------------------------
graph_embedding_length = 600
fingerprint_length = 1024

# MLP MODEL
class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, output_dim):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim1)
        self.fc2 = nn.Linear(hidden_dim1, hidden_dim2)
        self.fc3 = nn.Linear(hidden_dim2, output_dim)
        self.dropout = nn.Dropout(0.6)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        x = torch.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.fc3(x)
        return x
    
# GCN MODEL
class GCN(torch.nn.Module):
    def __init__(self, num_node_features, hidden_dim1, hidden_dim2, output_dim, chemberta):
        super(GCN, self).__init__()
        self.conv1 = GCNConv(num_node_features, hidden_dim1)
        self.conv2 = GCNConv(hidden_dim1, hidden_dim2)
        if chemberta:
            self.fc = torch.nn.Linear(hidden_dim2 + graph_embedding_length, output_dim)
        else:
            self.fc = torch.nn.Linear(hidden_dim2, output_dim)
        self.chemberta = chemberta

    def forward(self, x, edge_index, batch, edge_weight, graph_embedding):
        x = torch.relu(self.conv1(x.float(), edge_index, edge_weight))
        x = torch.relu(self.conv2(x, edge_index, edge_weight))
        x = global_add_pool(x, batch)
        if self.chemberta:
            x = torch.concat([x, torch.tensor(np.array(graph_embedding), dtype=torch.float).squeeze(1)], dim=1)
        x = self.fc(x)
        x = F.dropout(x, p=0.6, training=self.training)
        return x
 
# GIN MODEL 
class GIN(torch.nn.Module):
    def __init__(self, num_layers, hidden_dim, num_features, output_dim, chemberta):
        super(GIN, self).__init__()
        self.conv1 = GINConv(Sequential(
            Linear(num_features, hidden_dim),
            ReLU(),
            Linear(hidden_dim, hidden_dim),
            ReLU(),
            BN(hidden_dim)), train_eps=True)
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers - 1):
            self.convs.append(
                GINConv(Sequential(
                    Linear(hidden_dim, hidden_dim),
                    ReLU(),
                    Linear(hidden_dim, hidden_dim),
                    ReLU(),
                    BN(hidden_dim),
                ),
                        train_eps=True))
        if chemberta:
            self.lin1 = Linear(hidden_dim + graph_embedding_length, hidden_dim)
        else:
            self.lin1 = Linear(hidden_dim, hidden_dim)
        self.lin2 = Linear(hidden_dim, output_dim)
        self.chemberta = chemberta

    def forward(self, x, edge_index, batch, edge_weight, graph_embedding):
        x = self.conv1(x.float(), edge_index)
        for conv in self.convs:
            x = conv(x, edge_index)
        x = global_add_pool(x, batch)
        if self.chemberta:
            x = torch.concat([x, torch.tensor(np.array(graph_embedding), dtype=torch.float).squeeze(1)], dim=1)
        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.6, training=self.training)
        x = self.lin2(x)
        return x
    
# GAT MODEL
class GAT(torch.nn.Module):
    def __init__(self, num_features, hidden_channels, heads, output_dim, chemberta):
        super().__init__()
        torch.manual_seed(1234567)
        self.conv1 = GATConv(num_features, hidden_channels,heads)
        self.conv2 = GATConv(heads*hidden_channels, hidden_channels,heads)
        if chemberta:
            self.fc = torch.nn.Linear(heads*hidden_channels + graph_embedding_length, output_dim) 
        else:
            self.fc = torch.nn.Linear(heads*hidden_channels, output_dim)
        self.chemberta = chemberta

    def forward(self, x, edge_index, batch, edge_weight, graph_embedding):
        x = F.dropout(x, p=0.6, training=self.training)
        x = self.conv1(x.float(), edge_index, edge_weight)
        x = F.elu(x)
        x = F.dropout(x, p=0.6, training=self.training)

        x = self.conv2(x, edge_index, edge_weight)
        x = global_add_pool(x, batch)
        if self.chemberta:
            x = torch.concat([x, torch.tensor(np.array(graph_embedding), dtype=torch.float).squeeze(1)], dim=1)
        x = self.fc(x)

        return x