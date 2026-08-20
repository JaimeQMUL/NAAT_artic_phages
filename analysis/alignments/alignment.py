from src.biotools.msa import *


# Run MSA on references
files=['7Z3M.fasta', '9GBG.fasta', 'K35G_E198R.fasta', 'P04529.fasta']
RunMSA('1-DatabaseCuration/gold_standards', 'full_protein_sequence_alignment', files)