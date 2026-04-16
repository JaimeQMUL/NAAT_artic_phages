import pandas as pd
import subprocess
from src.biotools.fasta_tools import ExtractSequence
from src.helpers.data_wrangling import *
from src.visualisations.phylogenetic_trees import *
import os

# Align the filtered_database_fasta
os.chdir('uvsx/data')
subprocess.run('/opt/anaconda3/envs/cold/bin/famsa curated_database/filtered_curated_database.fasta  alignments/filtered_curated_database_alignment.fasta', shell=True)
print('Alignment completed successfully')

#Trim alignments from highly gapped regions
os.chdir('alignments')
subprocess.run('/opt/anaconda3/envs/cold/bin/trimal -automated1 -in filtered_curated_database_alignment.fasta -out filtered_curated_database_alignment_trimmed.fasta', shell=True)

#Get pairwise distance matrix of alignment from uvsx_t4
subprocess.run('/opt/anaconda3/envs/cold/bin/famsa -dist_export -square_matrix filtered_curated_database_alignment.fasta filtered_database_alignment_distance_matrix.csv', shell=True )

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
