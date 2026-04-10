import pandas as pd

from src.biotools.fasta_tools import ExtractSequence
from src.helpers.data_wrangling import *
from src.visualisations.phylogenetic_trees import *

#
# # # Cleaning short and long seqs from database
# # FilterFasta('uvsx/data/curated_database/cleaned_curated_database.fasta',
# #             'uvsx/data/curated_database/cleaned_metadata.csv',
# #             'uvsx/data/curated_database/filtered_curated_database.fasta')
#
# # Align the filtered_database_fasta
# # mafft filtered_curated_database.fasta > ../alignments/filtered_curated_database_alignment.fasta
#
# # Create newick tree
# # cd../alignments
# # FastTree filtered_curated_database_alignment.fasta > tree.nwk
#
# #Clean newick tree
# CleanNewickTree('uvsx/data/alignments/tree.nwk')
#
#
result_files=['uvsx/data/hmms/results/full_sequence_hmm_results.tblout',
         'uvsx/data/hmms/results/PF00154_custom_hmm_results.tblout',
         'uvsx/data/hmms/results/PF00154_hmm_results.tblout',
         'uvsx/data/hmms/results/PF21134_custom_hmm_results.tblout',
         'uvsx/data/hmms/results/PF21134_hmm_results.tblout']
#
mapping=MapAccessionsToHMMResults(result_files)
#
# PlotPhyloGroups('uvsx/data/alignments/tree_clean.nwk', mapping, 'uvsx/plots/trees')

PlotWholePhylogeny('uvsx/data/alignments/tree_clean.nwk','uvsx/plots/trees/whole_phylogeny.png')

PlotPhyloTreeGrouped('uvsx/data/alignments/tree_clean.nwk', mapping, 'uvsx/plots/trees/grouped_phylogeny.png')
