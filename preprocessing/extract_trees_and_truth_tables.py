from lxml import etree
import csv
import pandas as pd
from functools import reduce
import numpy as np
import sys

def descendant_loop(child, parent_name, level, df):
    """Recursive function to retrieve each descendant's parent name.

    Args:
        child: path to current node.
        parent_name: current node name
        level: depth level in the ontology tree (here between 1 and 6)
        df: dataframe where updates are made
    """
    if level == max_level+1: # No further descendants
        child.clear()
        return
    else:
        for c in child.iterfind('{http://www.hmdb.ca}descendants/{http://www.hmdb.ca}descendant'): # For each descendant of the current node
            c_name = c.find('{http://www.hmdb.ca}term').text # Retrieve the descendant's name
            df.at[c_name, 'Parent'] = parent_name # Update the relation between the current node and the descendant in the dataframe
            descendant_loop(c, c_name, level + 1, df) # Go one level deeper
            c.clear()

def ontology_loop(df):
    """Function to iterate through each metabolite and each one of its ontology nodes

    Args:
        df: dataframe where updates are made
    """
    context = etree.iterparse(xml_file, tag='{http://www.hmdb.ca}metabolite')
    for event, elem in context:
        for c0 in elem.iterfind('{http://www.hmdb.ca}ontology'): # Iterate through each metabolite
            for c1 in c0.iterfind('{http://www.hmdb.ca}root'): # Iterate through each root term
                c1_name = c1.find('{http://www.hmdb.ca}term').text
                df.at[c1_name, 'Parent'] = 'Root'
                descendant_loop(c1, c1_name, 2, df) # Recursive function to visit every node
        elem.clear()
        etree.cleanup_namespaces(elem)
    # Also eliminate now-empty references from the root node to elem
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    del context

def hmdbextract_ontology(in_file, ontology_terms, out_file):
    """Function to extract the ontology terms of each metabolite. This code has been modified from http://www.metabolomics-forum.com/index.php?topic=1588.0 

    Args:
        in_file: file from which the infos will be retrieved.
        ontology_terms: all the possible ontology terms in in_file
        out_file: file where the infos will be saved
    Return: 
        all_terms: nested dictionaries with infos for each level of each parent node and ID
    """
    ns = {'hmdb': 'http://www.hmdb.ca'}
    context = etree.iterparse(in_file, tag='{http://www.hmdb.ca}metabolite')
    csvfile = open(out_file, 'w')
    fieldnames = ['accession'] + ontology_terms # All the columns names

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    all_terms = {key: {} for key in np.arange(1, max_level+1)} # Initialize empty nested dictionaries
    for event, elem in context:
        accession = elem.xpath('hmdb:accession/text()', namespaces=ns)[0]
        # Level 1 has a slightly different skeleton. It needs to be parsed separately 
        try:
            level1_terms = elem.xpath('hmdb:ontology/hmdb:root/hmdb:term/text()', namespaces=ns)
            level1_ids= [-1 for i in level1_terms] 
        except:
            level1_terms = 'NA'
            level1_ids = 'NA'
        curr_dict = dict(zip(level1_terms, level1_ids))
        all_terms[1].update(curr_dict)
        all_parent_ids = level1_ids
        terms = level1_terms

        path = 'hmdb:ontology/hmdb:root/'
        for level in range(2,max_level+1): # For each level, retrieve the infos (term and parent ID)
            path = path + 'hmdb:descendants/hmdb:descendant/'
            try:
                curr_terms = elem.xpath(path + 'hmdb:term/text()', namespaces=ns)
                parent_ids = [int(i) for i in elem.xpath(path + 'hmdb:parent_id/text()', namespaces=ns)]
                all_parent_ids += parent_ids
                terms += curr_terms
            except:
                terms = 'NA'
                parent_ids = 'NA'
            curr_dict = dict(zip(curr_terms, parent_ids))
            all_terms[level].update(curr_dict)
        dict_to_write = {'accession': accession}
        for i in ontology_terms:
            dict_to_write.update({i: True if i in terms else False}) # Update the truth table
        writer.writerow(dict_to_write)
        # It's safe to call clear() here because no descendants will be
        # accessed
        elem.clear()
    # Also eliminate now-empty references from the root node to elem
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    del context
    return all_terms;

max_level = 6 # max level of descendant depth in the xml file 
xml_file = sys.argv[1] # path to xml file containing the metabolites infos
xml_idx = xml_file[38:40] # needs to be adapted based on the exact name of your files
truth_tables_path = sys.argv[2] # path to the folder where all the truth tables will be saved
trees_path = sys.argv[3] # path to the folder where all the trees will be saved

# Get the terms for all ontology nodes
print("Started to collect all the ontology nodes terms")
ontology_terms = pd.DataFrame()
ontology_loop(ontology_terms)
ontology_terms = ontology_terms.index.tolist()

df_tree = pd.DataFrame()
all_curr_terms = hmdbextract_ontology(xml_file, ontology_terms, truth_tables_path + "/truth_table{}.csv".format(xml_idx))
df_curr_tree = pd.DataFrame() # Start building the ontology tree
for key, value in all_curr_terms.items(): # Iterate through each level
    curr_df = pd.DataFrame.from_dict(value, orient='index', columns=['Parent ID']) # Transform each level dict into a dataframe
    curr_df['Level'] = key
    curr_df['ID'] = 0
    curr_df['Parent'] = ''
    curr_df = curr_df[['ID', 'Level', 'Parent ID', 'Parent']]
    df_tree = pd.concat([df_tree, curr_df]) # Concatenate the new dataframe

ontology_loop(df_tree)

if df_tree.size != 0:
    df_tree.at['Root', 'Parent ID'] = -2 
    df_tree.at['Root', 'Level'] = 0
    df_tree['Parent ID'] = df_tree['Parent ID'].astype(int)
    df_tree['Level'] = df_tree['Level'].astype(int)
    for i in range(6,0,-1): # Iterate through each level
        result = df_tree[df_tree['Level'] == i]
        # For each term, update the Parent ID, the Parent, and the ID
        for i in result.iterrows():
            parent_id = df_tree.loc[i[0], 'Parent ID']
            parent = df_tree.loc[i[0], 'Parent']
            df_tree.at[parent, 'ID'] = parent_id

    next_available_id = df_tree['ID'].sort_values(ascending=False).iloc[0] + 1 # Find next available ID for term without an ID (level6)
    for idx, row in df_tree[df_tree['ID'] == 0].iterrows():
        df_tree.at[idx, 'ID'] = next_available_id
        next_available_id += 1
    df_tree = df_tree.sort_values('ID')

    df_tree['ID'] = df_tree['ID'].astype(int)
    df_tree = df_tree.reset_index()
    df_tree.columns.values[0] = 'Node'
    df_tree.to_csv(trees_path + "/tree{}.csv".format(xml_idx))
    print("Saved the ontology tree to tree{}.csv".format(xml_idx))
else:
    print("No ontology for this set of metabolites")