from rdkit import Chem
import numpy as np
import networkx as nx
from torch_geometric.data import Data
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler
import pandas as pd
import torch
import sys

class GraphDataset(Dataset):
    def __init__(self, data_list):
        self.dataset = data_list
    
    def __len__(self):
        return len(self.dataset)
    
    def __getitem__(self, idx):
        return self.dataset[idx]

def mol_to_graph(mol):
    """Convert RDKit molecule to networkx graph

    Args:
        mol: RDKit molecule
    Return: 
        G: networkx graph
    """
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), atomic_num=atom.GetAtomicNum())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondTypeAsDouble())
    return G

def graph_to_tensors(G):
    """Convert networkx graph to PyTorch tensors

    Args:
        G: networkx graph
    Return: 
        node_features: tensor containing information on nodes (3D coordinates and atomic numbers) of G
        edge_features: tensor containing information on bond types of G
        adjacency matrix: tensor containing adjacency information of G
    """
    node_features = np.array([G.nodes[i]['atomic_num'] for i in G.nodes()], dtype=np.float32)
    edge_features = np.array([G.edges[edge]['bond_type'] for edge in G.edges()], dtype=np.float32)
    adjacency_matrix = nx.to_numpy_array(G)
    return torch.FloatTensor(node_features), torch.FloatTensor(edge_features), torch.FloatTensor(adjacency_matrix)

sdf_file = sys.argv[1] # Path to sdf file containing structures
output_file = sys.argv[2] # Path to outputs file
saving_folder = sys.argv[3] # Path to folder where files will be saved
selected_ids = np.load(sys.argv[4], allow_pickle=True) # IDs to keep (usually only quantified metabolites)
selected_nodes = np.concatenate((['accession'], np.load(sys.argv[5], allow_pickle=True))) # Nodes to keep (given a level 1 node and std >0.1)
extra_info = pd.read_csv(sys.argv[6]) # Extra information to use as input such as average molecular weight
extra_info_ids = extra_info['accession'].values
mol_supplier = Chem.SDMolSupplier(sdf_file)

num_graphs = len(mol_supplier)
max_num_nodes = 0  # Track the maximum number of nodes among all graphs
max_num_edges = 0  # Track the maximum number of edges among all graphs

# Initialize arrays to store information for all graphs
node_features_array = []
edge_features_array = []
adjacency_matrix_array = []
valid_ids = []
num_of_2d_mols = 0
possible_atoms = set()

for i, mol in enumerate(mol_supplier):
    if i%10000==0:
        print("Processed ", i, " Molecules")
    if mol is not None:
        curr_id = mol.GetProp('DATABASE_ID')
        # Check that the metabolite is quantified and we have extra info about it
        if curr_id in selected_ids and curr_id in extra_info_ids:
            # Get nodes and edges of molecule
            conf = mol.GetConformer()
            num_nodes = mol.GetNumAtoms()
            node_features = torch.empty((num_nodes,3))
            for atom_idx in range(num_nodes):
                atom = mol.GetAtomWithIdx(atom_idx)
                possible_atoms.add(atom.GetAtomicNum())
                x, y, z = conf.GetAtomPosition(atom_idx)
                node_features[atom_idx] = torch.tensor([x, y, z])
            # Count number of 2D structures
            if node_features[:,2].sum()==0:
                num_of_2d_mols += 1

            # Standardize 3D coordinates
            scaler = StandardScaler()
            scaler.fit(node_features)
            node_features = scaler.transform(node_features)

            graph = mol_to_graph(mol) # Convert molecule into a graph
            node_atomic_num, edge_features, adjacency_matrix = graph_to_tensors(graph)
        
            node_features = torch.cat((node_atomic_num.unsqueeze(0).t(), node_features), dim=1)
            num_edges = len(edge_features)

            # Remove metabolites with only one edge
            if num_edges > 1:
                valid_ids.append(curr_id)
                max_num_nodes = max(max_num_nodes, num_nodes) # Update maximum number of nodes 
                max_num_edges = max(max_num_edges, num_edges) # Update maximum number of edges 
            
                node_features_array.append(node_features)
                edge_features_array.append(edge_features)
                adjacency_matrix_array.append(adjacency_matrix)

possible_atoms = torch.tensor(np.sort(list(possible_atoms))) # List of possible atoms (for one-hot encoding)
np.savetxt("used_atoms.txt", np.sort(np.array(list(possible_atoms), dtype=int)))
num_graphs = len(valid_ids) # Keep only valid graphs
print("Number of graphs:", num_graphs)

node_features_array_padded = np.zeros((num_graphs, max_num_nodes, 3+len(possible_atoms)), dtype=np.float16)
edge_features_array_padded = np.zeros((num_graphs, max_num_edges), dtype=np.float16)
adjacency_features = np.zeros((num_graphs, 2, max_num_edges), dtype=np.int16)
all_ys = np.zeros((num_graphs, len(selected_nodes)-1))
all_graphs = []

output = pd.read_csv(output_file)
output = output[selected_nodes] # Filter to keep only selected nodes (ontology terms)
filtered_output = output[output['accession'].isin(valid_ids)] # Filter to keep only valid metabolites

print("Number of molecules without 3D coordinates: ", num_of_2d_mols)

for i in range(num_graphs):
    if i%1000==0:
        print(f"Processed {i}/{num_graphs} Graphs")
    # Get graph elements to create Data object
    curr_node_features_array = node_features_array[i]
    curr_edge_features_array = edge_features_array[i]
    curr_adjacency_features = adjacency_matrix_array[i]
    
    num_nodes = len(curr_node_features_array)
    num_edges = len(curr_edge_features_array)

    edge_features_array_padded[i, :num_edges] = curr_edge_features_array

    row, col = curr_adjacency_features.triu().nonzero(as_tuple=True) # Convert adjacency matrix into smaller array 
    updated_adj = torch.stack([row, col], dim=0)
    adjacency_features[i, :, :num_edges] = updated_adj

    curr_y = filtered_output.loc[filtered_output['accession']==valid_ids[i]].values[0,1:]
    all_ys[i,:] = curr_y
    curr_vector = curr_node_features_array[:,0]

    one_hot_encoded_matrix = (curr_vector.unsqueeze(1) == possible_atoms).float() # One hot encode atomic number
    curr_extra_info = extra_info[extra_info['accession']==valid_ids[i]]['Avg mol. weight'].values
    curr_node_features_array = torch.cat([curr_node_features_array[:,1:], one_hot_encoded_matrix], dim=1)
    node_features_array_padded[i, :num_nodes,:] = curr_node_features_array

    # Add new graph to list of all graphs
    all_graphs.append(Data(x= curr_node_features_array, edge_index=updated_adj, edge_attr=curr_edge_features_array,
                           y=np.nan_to_num(curr_y.astype(bool)), extra_info = curr_extra_info))

dataset = GraphDataset(all_graphs) # Convert  list of graphs into PyTorch Geometric object

# Save all processed files
print("Node features array shape:", node_features_array_padded.shape)
np.save(saving_folder + "/node_features.npy", node_features_array_padded) # Save node features (3D coordinates + one-hot encoded atomic number) of each graph
print("Edge features array shape:", edge_features_array_padded.shape)
np.save(saving_folder + "/edge_features.npy", edge_features_array_padded) # Save edge features (bond type) of each graph
print("Adjacency features array shape:", adjacency_features.shape)
np.save(saving_folder + "/adjacency_features.npy", adjacency_features) # Save adjacency features (nodes of bond) of eah graph
filtered_output.to_csv(saving_folder + "/filtered_outputs.csv", index=False) # Save the filtered outputs to keep only valid IDs
torch.save(dataset, saving_folder + "/dataset.pt") # Save dataset as a PyTorch Geometric object
print("Number of valid ids: ", len(valid_ids)) 
np.save(saving_folder + "/valid_ids", valid_ids) # Save ids of all graphs