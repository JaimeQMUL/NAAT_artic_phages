from src.helpers.data_wrangling import *

# Filtering short and long seqs from database
filtered_accs=FilterFasta('uvsx/data/curated_database/cleaned_curated_database.fasta',
            'uvsx/data/curated_database/cleaned_metadata.csv',
            'uvsx/data/curated_database/filtered_curated_database.fasta')


# Filtering metadata to only have seqs in filtered data
#  create new csv with only seqs from filtered fasta
print(len(filtered_accs))