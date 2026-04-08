
################################################################################################################################################
# Cleaning short seqs from database
################################################################################################################################################
import pandas as pd
import numpy as np
from src.biotools.fasta_tools import ExtractSequence
from src.helpers.data_wrangling import *


df=pd.read_csv('uvsx/data/curated_database/cleaned_metadata.csv', low_memory=True)
lengths=[]
filtered=[]
for acc in df['accession']:
    seq=ExtractSequence(acc,'uvsx/data/curated_database/cleaned_curated_database.fasta' )
    length=len(seq)
    lengths.append(length)

mean_length = np.mean(lengths)
std_length = np.std(lengths)

lower = mean_length - 2*std_length
upper = mean_length + 2*std_length

for acc in df['accession']:
    seq=ExtractSequence(acc,'uvsx/data/curated_database/cleaned_curated_database.fasta')
    length=len(seq)
    if length > lower and length < upper:
        filtered.append(acc)




WriteCleanedFasta(filtered, 'uvsx/data/curated_database/filtered_curated_database.fasta', 'uvsx')

# Align the filtered_database_fasta
# mafft filtered_curated_database.fasta > ../alignments/filtered_curated_database_alignment.fasta