import pandas as pd
import subprocess
from src.biotools.fasta_tools import ExtractSequence
from src.helpers.data_wrangling import *
from src.visualisations.phylogenetic_trees import *
import os


# ================================================== FILTERED DATABASE =================================================
# # File for analysis
# file='filtered_curated_database'
#
# # Create neccessary directories
# os.chdir('uvsx/data')
# os.makedirs(f'alignments/{file}', exist_ok=True)
# os.makedirs(f'iqtree/{file}', exist_ok=True)
#
# # Align the database fasta
# subprocess.run(f'/opt/anaconda3/envs/cold/bin/famsa curated_database/{file}.fasta  alignments/{file}/{file}.aln.fasta', shell=True)
# print('Alignment completed successfully')
#
# #Trim alignments from highly gapped regions
# os.chdir(f'alignments/{file}')
# subprocess.run(f'/opt/anaconda3/envs/cold/bin/trimal -automated1 -in {file}.aln.fasta -out {file}.aln.trm.fasta', shell=True)
#
# #Get pairwise distance matrix of alignment from uvsx_t4
# subprocess.run(f'/opt/anaconda3/envs/cold/bin/famsa -dist_export -square_matrix {file}.aln.trm.fasta {file}_distance_matrix.csv', shell=True )
#
# #Get Phylogenetic tree from iqtree
# os.chdir('../..')
# subprocess.run(f'/opt/anaconda3/envs/cold/bin/iqtree3 -s alignments/{file}/{file}.aln.trm.fasta --prefix iqtree/{file}/{file}_tree', shell=True)

# ======================================================================================================================
# ====================================================== TOP HITS ======================================================
# File for analysis
file='top_hits'

# Create neccessary directories
os.chdir('uvsx/data')
os.makedirs(f'alignments/{file}', exist_ok=True)
os.makedirs(f'iqtree/{file}', exist_ok=True)

# Align the database fasta
subprocess.run(f'/opt/anaconda3/envs/cold/bin/famsa curated_database/{file}.fasta  alignments/{file}/{file}.aln.fasta', shell=True)
print('Alignment completed successfully')

#Trim alignments from highly gapped regions
os.chdir(f'alignments/{file}')
subprocess.run(f'/opt/anaconda3/envs/cold/bin/trimal -automated1 -in {file}.aln.fasta -out {file}.aln.trm.fasta', shell=True)

#Get pairwise distance matrix of alignment from uvsx_t4
subprocess.run(f'/opt/anaconda3/envs/cold/bin/famsa -dist_export -square_matrix {file}.aln.trm.fasta {file}_distance_matrix.csv', shell=True )

#Get Phylogenetic tree from iqtree
os.chdir('../..')
subprocess.run(f'/opt/anaconda3/envs/cold/bin/iqtree3 -s alignments/{file}/{file}.aln.trm.fasta --prefix iqtree/{file}/{file}_tree', shell=True)







# # Create newick tree
# subprocess.run(f'FastTree {align_dir}filtered_curated_database_alignment.fasta > {align_dir}tree.nwk', shell=True)
#
# #Clean newick tree
# CleanNewickTree('uvsx/data/alignments/tree.nwk')
#
#
# result_files=['uvsx/data/hmms/results/full_sequence_hmm_results.tblout',
#          'uvsx/data/hmms/results/PF00154_custom_hmm_results.tblout',
#          'uvsx/data/hmms/results/PF00154_hmm_results.tblout',
#          'uvsx/data/hmms/results/PF21134_custom_hmm_results.tblout',
#          'uvsx/data/hmms/results/PF21134_hmm_results.tblout']
#
# mapping=MapAccessionsToHMMResults(result_files)
#
# PlotPhyloGroups('uvsx/data/alignments/tree_clean.nwk', mapping, 'uvsx/plots/trees')
#
# PlotWholePhylogeny('uvsx/data/alignments/tree_clean.nwk','uvsx/plots/trees/whole_phylogeny.png')
#
# PlotPhyloTreeGrouped('uvsx/data/alignments/tree_clean.nwk', mapping, 'uvsx/plots/trees/grouped_phylogeny.png')
