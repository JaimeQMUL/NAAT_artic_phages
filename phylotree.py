import pandas as pd
import subprocess
from src.biotools.fasta_tools import ExtractSequence
from src.helpers.data_wrangling import *
from src.visualisations.phylogenetic_trees import *

#Establish directories
ref_dir='uvsx/data/references/'
align_dir='uvsx/data/alignments/'
data_dir='uvsx/data/curated_database/'


# Align the filtered_database_fasta
# Add references to the database
subprocess.run(f'cat {ref_dir}P04529.fasta {ref_dir}7Z3M.fasta {ref_dir}9GBG.fasta {ref_dir}K35G_E198R.fasta >> {data_dir}filtered_curated_database.fasta', shell=True)


# Run alignment

subprocess.run(f'/opt/anaconda3/envs/cold/bin/mafft --localpair --maxiterate 1000 {data_dir}filtered_curated_database.fasta > {align_dir}filtered_curated_database_alignment.fasta', shell=True)
print('Alignment completed successfully')

#Trim alignments from highly gapped regions


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
