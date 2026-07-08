from src.helpers.data_wrangling import *
import os, subprocess

# Filtering short and long seqs from database
filtered_accs=FilterFasta('uvsx/data/curated_database/cleaned_curated_database.fasta',
            'uvsx/data/curated_database/cleaned_metadata.csv',
            'uvsx/data/curated_database/filtered_curated_database.fasta')

# add Gold standards to filtered database
os.chdir('uvsx/data/references')
subprocess.run(f'cat 7Z3M.fasta 9GBG.fasta K35G_E198R.fasta >> ../curated_database/filtered_curated_database.fasta', shell=True)



# Filtering metadata to only have seqs in filtered data
#  create new csv with only seqs from filtered fasta

df = pd.DataFrame(filtered_accs, columns=["accession"])
df.to_csv("uvsx/data/curated_database/filtered_metadata.csv", index=False)

print(len(filtered_accs))