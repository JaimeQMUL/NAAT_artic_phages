# Clean data removing all redundant sequences.
# These are sequences that have been completely duplicated, same accession same sequence.
# As well as removing sequences that have different accessions, but the same sequence and are from the same organism

#Imports
from src.helpers.data_wrangling import *
import pandas as pd

# Get all the unique accessions from the initial database search
unique_accessions=GetUniqueAccessions('uvsx')

# Find the accessions of non-redudant proteins in the initial database
non_dupes=RemoveRedundancies(unique_accessions, 'uvsx')

#Write to csv for later
df = pd.DataFrame(non_dupes, columns=['accession'])
df.to_csv('uvsx/data/curated_database/cleaned_metadata.csv', index=False)

#Adding p04529 to non dupes. Replacing the accession being used instead
test='P04529' in non_dupes
print(f'P04529 in non_dupes: {test}')
try:
    i = non_dupes.index('ADJ39758')
    non_dupes[i] = 'P04529'
except ValueError:
    pass

WriteCleanedFasta(non_dupes, 'uvsx')

#Check there are no discrepencies
discrepencies=0
dis_list=[]
for acc in non_dupes:
    seq=ExtractSequence(acc, 'uvsx/data/curated_database/cleaned_curated_database.fasta')
    test_seq=ExtractSequence(acc, 'uvsx/data/curated_database/uniprot_query_search.fasta')
    if test_seq=='':
        test_seq=ExtractSequence(acc, 'uvsx/data/curated_database/ncbi_query_search.fasta')
    if seq != test_seq:
        discrepencies+=1
        dis_list.append(acc)
print(f'discrepencies: {discrepencies}')



