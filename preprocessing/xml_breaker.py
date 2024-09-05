from lxml import etree
import csv
import pandas as pd
from functools import reduce
import sys
import os

def hmdbextract(in_file, out_file):
  """Function to extract all the information (except ontology terms) of each metabolite. This code has been modified from http://www.metabolomics-forum.com/index.php?topic=1588.0 

    Args:
        in_file: file from which the infos will be retrieved.
        out_file: file where the infos will be saved
    """
  ns = {'hmdb': 'http://www.hmdb.ca'}
  context = etree.iterparse(in_file, tag='{http://www.hmdb.ca}metabolite')
  csvfile = open(out_file, 'w')
  fieldnames = ['accession', 'status',
                'Avg mol. weight', 'iupac_name', 'chemical_formula', 'InChIKey', 'InChI', 'smiles', 
                'kingdom',  'direct_parent', 'super_class', 'class', 'sub_class', 'molecular_framework', 'alternative_parents', 'substituents', 'external_descriptors',
                'state', 'exp_melting_point','exp_boiling_point', 'exp_water_solubility', 'exp_logP',
                'cellular_locations', 'biospecimen_locations', 'tissue_locations', 'pathway_name', 'pathway_smpdb_id', 'pathway_kegg_map_id',
                'drugbank','chebi_id', 'pubchem', 'phenol_explorer_compound_id','food','knapsack', 'chemspider', 'kegg', 'biocyc_id','bigg','metlin_id','pdb_id', 'vmh_id', 'fbonto_id',
                'diseases', 'diseases_omim', 'enzymes']
  writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
  writer.writeheader()
  
  for event, elem in context:
    # RECORD INFORMATION
    accession = elem.xpath('hmdb:accession/text()', namespaces=ns)[0]
    try:
        status = elem.xpath('hmdb:status/text()', namespaces=ns)[0].encode('utf-8')
    except:
        status = 'NA'
    # METABOLITE IDENTIFICATION
    try:
        average_molecular_weight = elem.xpath('hmdb:average_molecular_weight/text()', namespaces=ns)[0]
    except:
        average_molecular_weight = 'NA'
    try:
        iupac_name = elem.xpath('hmdb:iupac_name/text()', namespaces=ns)[0].encode('utf-8')
    except:
        iupac_name = 'NA'
    try:
        chemical_formula = elem.xpath('hmdb:chemical_formula/text()', namespaces=ns)[0]
    except:
        chemical_formula = 'NA'
    try:
        inchikey = elem.xpath('hmdb:inchikey/text()', namespaces=ns)[0]
    except:
        inchikey = 'NA'
    try:
        inchi = elem.xpath('hmdb:inchi/text()', namespaces=ns)[0]
    except:
        inchi = 'NA'
    try:
        smiles = elem.xpath('hmdb:smiles/text()', namespaces=ns)[0]
    except:
        smiles = 'NA'
    
    # TAXONOMY
    try:
        kingdom = elem.xpath('hmdb:taxonomy/hmdb:kingdom/text()', namespaces=ns)[0]
    except:
        kingdom = 'NA'
    try:
        direct_parent = elem.xpath('hmdb:taxonomy/hmdb:direct_parent/text()', namespaces=ns)[0]
    except:
        direct_parent = 'NA'
    try:
        super_class = elem.xpath('hmdb:taxonomy/hmdb:super_class/text()', namespaces=ns)[0]
    except:
        super_class = 'NA'
    try:
        classorg = elem.xpath('hmdb:taxonomy/hmdb:class/text()', namespaces=ns)[0]
    except:
        classorg = 'NA'
    try:
        sub_class = elem.xpath('hmdb:taxonomy/hmdb:sub_class/text()', namespaces=ns)[0]
    except:
        sub_class = 'NA'
    try:
        molecular_framework = elem.xpath('hmdb:taxonomy/hmdb:molecular_framework/text()', namespaces=ns)[0]
    except:
        molecular_framework = 'NA'
    try:
        alternative_parents = elem.xpath('hmdb:taxonomy/hmdb:alternative_parents/hmdb:alternative_parent/text()', namespaces=ns)
        if not alternative_parents:
            alternative_parents = 'NA'
    except:
        alternative_parents = 'NA'
    try:
        substituents = elem.xpath('hmdb:taxonomy/hmdb:substituents/hmdb:substituent/text()', namespaces=ns)
        if not substituents:
            substituents = 'NA'
    except:
        substituents = 'NA'
    try:
        external_descriptors = elem.xpath('hmdb:taxonomy/hmdb:external_descriptors/hmdb:external_descriptor/text()', namespaces=ns)
        if not external_descriptors:
            external_descriptors = 'NA'
    except:
        external_descriptors = 'NA'
      
    # PHYSICAL PROPERTIES
    try:
        state = elem.xpath('hmdb:state/text()', namespaces=ns)[0]
    except:
        state = 'NA'
    try:
        experimental_properties_kind = elem.xpath('hmdb:experimental_properties/hmdb:property/hmdb:kind/text()', namespaces=ns)
    except:
        experimental_properties_kind = 'NA'
    try:
        experimental_properties_values = elem.xpath('hmdb:experimental_properties/hmdb:property/hmdb:value/text()', namespaces=ns)
    except:
        experimental_properties_values = 'NA'

    exp_melting_point = 'NA'
    exp_boiling_point = 'NA'
    exp_water_solubility = 'NA'
    exp_logP = 'NA'
    for idx, k in enumerate(experimental_properties_kind):
        if k == 'melting_point':
            exp_melting_point = experimental_properties_values[idx]
        elif k == 'boiling_point':
            exp_boiling_point = experimental_properties_values[idx]
        elif k == 'water_solubility':
            exp_water_solubility = experimental_properties_values[idx]
        elif k == 'logp':
            exp_logP = experimental_properties_values[idx]
    
    # BIOLOGICAL PROPERTIES
    try:
        cellular_locations = elem.xpath('hmdb:biological_properties/hmdb:cellular_locations/hmdb:cellular/text()', namespaces=ns)
        if not cellular_locations:
            cellular_locations = 'NA'
    except:
        cellular_locations = 'NA'

    try:
        biospecimen_locations = elem.xpath('hmdb:biological_properties/hmdb:biospecimen_locations/hmdb:biospecimen/text()', namespaces=ns)
        if not biospecimen_locations:
            biospecimen_locations = 'NA'
    except:
        biospecimen_locations = 'NA'

    try:
        tissue_locations = elem.xpath('hmdb:biological_properties/hmdb:tissue_locations/hmdb:tissue/text()', namespaces=ns)
        if not tissue_locations:
            tissue_locations = 'NA'
    except:
        tissue_locations = 'NA'
    
    try:
        pathway_name = elem.xpath('hmdb:biological_properties/hmdb:pathways/hmdb:pathway/hmdb:name/text()', namespaces=ns)
        if not pathway_name:
            pathway_name = 'NA'
    except:
        pathway_name = 'NA'
    try:
        pathway_smpdb_id = elem.xpath('hmdb:biological_properties/hmdb:pathways/hmdb:pathway/hmdb:smpdb_id/text()', namespaces=ns)
        if not pathway_smpdb_id:
            pathway_smpdb_id = 'NA'
    except:
        pathway_smpdb_id = 'NA'
    try:
        pathway_kegg_map_id = elem.xpath('hmdb:biological_properties/hmdb:pathways/hmdb:pathway/hmdb:kegg_map_id/text()', namespaces=ns)
        if not pathway_kegg_map_id:
            pathway_kegg_map_id = 'NA'
    except:
        pathway_kegg_map_id = 'NA'

    # EXTERNAL LINKS
    try:
        drugbank = elem.xpath('hmdb:drugbank_id/text()', namespaces=ns)[0]
    except:
        drugbank = 'NA'
    try:
        chebi_id = elem.xpath('hmdb:chebi_id/text()', namespaces=ns)[0]
    except:
        chebi_id = 'NA'
    try:
        pubchem = elem.xpath('hmdb:pubchem_compound_id/text()', namespaces=ns)[0]
    except:
        pubchem = 'NA'
    try:
        phenol_explorer_compound_id = elem.xpath('hmdb:phenol_explorer_compound_id/text()', namespaces=ns)[0]
    except:
        phenol_explorer_compound_id = 'NA'
    try:
        food = elem.xpath('hmdb:foodb_id/text()', namespaces=ns)[0]
    except:
        food = 'NA'
    try:
        knapsack = elem.xpath('hmdb:knapsack_id/text()', namespaces=ns)[0]
    except:
        knapsack = 'NA'
    try:
        chemspider = elem.xpath('hmdb:chemspider_id/text()', namespaces=ns)[0]
    except:
        chemspider = 'NA'
    try:
        kegg = elem.xpath('hmdb:kegg_id/text()', namespaces=ns)[0]
    except:
        kegg = 'NA'
    try:
        biocyc_id = elem.xpath('hmdb:biocyc_id/text()', namespaces=ns)[0]
    except:
        biocyc_id = 'NA'
    try:
        bigg = elem.xpath('hmdb:bigg_id/text()', namespaces=ns)[0]
    except:
        bigg = 'NA'
    try:
        metlin_id = elem.xpath('hmdb:metlin_id/text()', namespaces=ns)[0]
    except:
        metlin_id = 'NA'
    try:
        pdb_id = elem.xpath('hmdb:pdb_id/text()', namespaces=ns)[0]
    except:
        pdb_id = 'NA'
    try:
        vmh_id = elem.xpath('hmdb:vmh_id/text()', namespaces=ns)[0]
    except:
        vmh_id = 'NA'
    try:
        fbonto_id = elem.xpath('hmdb:fbonto_id/text()', namespaces=ns)[0]
    except:
        fbonto_id = 'NA'
        
    # ASSOCIATED DISORDERS AND DISEASES
    try:
        diseases = elem.xpath('hmdb:diseases/hmdb:disease/hmdb:name/text()', namespaces=ns)
        if not diseases:
            diseases = 'NA'
    except:
        diseases = 'NA'
    try:
        diseases_omim = elem.xpath('hmdb:diseases/hmdb:disease/hmdb:omim_id/text()', namespaces=ns)
        if not diseases_omim:
            diseases_omim = 'NA'
    except:
        diseases_omim = 'NA'

    # ENZYMES
    try:
        enzymes = elem.xpath('hmdb:protein_associations/hmdb:protein/hmdb:protein_accession/text()', namespaces=ns)
        if not enzymes:
            enzymes = 'NA'
    except:
        enzymes = 'NA'
    
    writer.writerow({'accession': accession, 'status': status, 
                     'Avg mol. weight': average_molecular_weight, 'iupac_name': iupac_name, 'chemical_formula': chemical_formula, 'InChIKey': inchikey, 'InChI': inchi, 'smiles': smiles,
                     'kingdom': kingdom, 'direct_parent': direct_parent, 'super_class': super_class, 'class': classorg, 'sub_class': sub_class, 'molecular_framework': molecular_framework, 'alternative_parents': alternative_parents, 'substituents': substituents, 'external_descriptors': external_descriptors,
                     'state': state, 'exp_melting_point': exp_melting_point,'exp_boiling_point': exp_boiling_point, 'exp_water_solubility': exp_water_solubility, 'exp_logP': exp_logP,
                     'cellular_locations': cellular_locations, 'biospecimen_locations': biospecimen_locations, 'tissue_locations': tissue_locations, 'pathway_name': pathway_name, 'pathway_smpdb_id': pathway_smpdb_id, 'pathway_kegg_map_id': pathway_kegg_map_id,
                     'drugbank': drugbank,'chebi_id': chebi_id,'pubchem': pubchem,'phenol_explorer_compound_id':phenol_explorer_compound_id, 'food': food,'knapsack': knapsack, 'chemspider': chemspider,'kegg': kegg, 'biocyc_id': biocyc_id, 'bigg':bigg, 'metlin_id': metlin_id, 'pdb_id':pdb_id, 'vmh_id': vmh_id, 'fbonto_id': fbonto_id,
                     'diseases': diseases, 'diseases_omim': diseases_omim, 'enzymes': enzymes})
    # It's safe to call clear() here because no descendants will be
    # accessed
    elem.clear()
# Also eliminate now-empty references from the root node to elem
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]
  del context
  return;

xml_file = sys.argv[1] # path to xml file containing the metabolites infos
truth_tables_path = sys.argv[2] # path to the folder where all the truth tables will be saved
trees_path = sys.argv[3] # path to the folder where all the trees will be saved

print("Started to merge the ontology trees")
# Update the ontology tree
df_tree = pd.DataFrame(columns=['Node'])

for filename in os.listdir(trees_path):
    curr_df = pd.read_csv(trees_path + '/' + filename, index_col=[0])
    df_tree = pd.concat([df_tree, curr_df]);
grouped = df_tree.groupby(by='Node')
df_tree = grouped.first().reset_index()
df_tree['Parent ID'] = df_tree['Parent ID'].astype(int)
df_tree['Level'] = df_tree['Level'].astype(int)
df_tree['ID'] = df_tree['ID'].astype(int)
df_tree.to_csv("ontology_tree.csv", index=False)
print("Saved the final ontology tree")

print("Started to merge the truth tables")
fieldnames = []
for filename in os.listdir(truth_tables_path):
    with open(truth_tables_path + '/' + filename, "r", newline="") as f_in:
        reader = csv.reader(f_in)
        headers = next(reader)
        for h in headers:
          if h not in fieldnames:
            fieldnames.append(h)
    
with open("ontology_truth_table.csv", "w", newline="") as f_out:
    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
    writer.writeheader() #this is the addition.       
    for filename in os.listdir(truth_tables_path):
        with open(truth_tables_path + '/' + filename, "r", newline="") as f_in:
            reader = csv.DictReader(f_in)  # Uses the field names in this file
            for line in reader:
                writer.writerow(line)
print("Saved the final truth table")

print("Started to collect infos about each metabolite")
hmdbextract(xml_file, 'metabolites_infos.csv') # Get all infos about each metabolite
print("Saved all the metabolite infos to metabolites_infos.csv")