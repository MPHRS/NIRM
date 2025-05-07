#Gene ID → GCA genome accession → Strain name
import pandas as pd
from Bio import Phylo
import re
import os
from io import StringIO

input_dir = "./trees_5_variab"#"./trees_5_shannon"
output_dir = "./modified_trees_variab" #"./modified_trees"
os.makedirs(output_dir, exist_ok=True)

# Load data_summary.tsv
data_summary = pd.read_csv('data_summary.tsv', sep='\t')
sample_to_strain = dict(zip(data_summary['Assembly Accession'], data_summary['Assembly Name']))

# Load cds_ortho.tsv
cds_ortho = pd.read_csv('cds_ortho.proteinortho.tsv', sep='\t')

# Map gene IDs to GCA IDs
gene_to_gca = {}
for _, row in cds_ortho.iterrows():
    for gca_col in cds_ortho.columns[3:]:
        genes = str(row[gca_col]).split(',')
        gca = gca_col.split('_cds_')[0]
        for gene in genes:
            gene_id = gene.strip()
            if gene_id not in ('*', 'nan', ''):
                gene_to_gca[gene_id] = gca

'''
# no iteration script
# Load original Newick tree
with open('aligned_1.fna_backup.treefile', 'r') as file:
    original_tree_str = file.read()
'''

# Replace gene IDs with strain names
def replace_gene_with_strain(match):
    gene_id = match.group(1)
    gca_id = gene_to_gca.get(gene_id, "Unknown")
    return sample_to_strain.get(gca_id, gca_id)

'''
# no iteration script
# Replace only leaf names, preserving branch lengths
modified_tree_str = re.sub(r'(lcl\|[^:(),]+)', replace_gene_with_strain, original_tree_str)

# Save the modified Newick tree
with open('strain_tree.nwk', 'w') as file:
    file.write(modified_tree_str)
'''

for file in os.listdir(input_dir):
    if file.endswith(".treefile"):
        input_path = os.path.join(input_dir, file)
        output_path = os.path.join(output_dir, file.replace(".treefile", "_strain.treefile"))
        
        # Load and process the original tree
        with open(input_path, 'r') as f:
            original_tree_str = f.read()
        
        # Replace gene IDs with strain names
        modified_tree = re.sub(r"(lcl\|[\w\.]+)", replace_gene_with_strain, original_tree_str)
        
        # Save modified tree
        with open(output_path, 'w') as f:
            f.write(modified_tree)
        
        print(f"Processed {file} → {os.path.basename(output_path)}")