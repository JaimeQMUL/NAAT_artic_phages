import pandas as pd
from src.biotools.fasta_tools import *
from src.helpers.data_wrangling import *
import os
# Initialise The scoring csv.
# Will contain the filtered accessions

hits=pd.read_csv('../1-DatabaseCuration/filtered/metadata/filtered_metadata.csv')
print(f'Number of accessions: {len(hits)}')

# Adding the Results from HMMs to hits
results_dir='../2-HiddenMarkovModels/results/'
results=os.listdir(results_dir)
mapping={}
for file in results:
    names=file.split('/')
    name=names[-1].split('.')[0]
    found=FindAccessionsInHMMResults(f'{results_dir}{file}')
    print(f'{name}: Number of hits {len(found)}')
    mapping[name]=found

for key in mapping.keys():
    binary = []

    for acc in hits['accession']:
        if acc in mapping[key]:
            score = 1
        else:
            score = 0
        binary.append(score)

    hits[key] = binary


# Adding results from Walker A and B regex search

# Finding walker A
a_results, a_present=MotifDatabaseSearch(hits['accession'], '../1-DatabaseCuration/clean/sequences/cleaned_curated_database.fasta', r"[GA]....[GF]K[TS]")
print(f'Walker A found {sum(a_present)}/{len(hits)}')
hits['Walker_A']=a_present

# Finding walker B
b_results, b_present=MotifDatabaseSearch(hits['accession'], '../1-DatabaseCuration/clean/sequences/cleaned_curated_database.fasta', r"[AVILMFWY]{4}DE?")
print(f'Walker B found {sum(b_present)}/{len(hits)}')
hits['Walker_B']=b_present

#Writing it to a csv for later
hits.to_csv('scoring_results/filtered_scoring.csv', index=False)

# Selecting accessions that hit for every check except the PF21134 hmms
target_cols = [
    'PF00154_results',
    'PF00154_domain_results',
    'gold_standards_results',
    'Walker_A',
    'Walker_B'
]

top_hits = hits[(hits[target_cols] == 1).all(axis=1)]['accession']

# Writing these sequences to a fasta
WriteCleanedFasta(top_hits, 'top_hits/sequences/top_hits.fasta')
print(f'{len(top_hits)}/{len(hits)} Written to fasta')

# Writing accessions to a csv for later
top_hits.to_csv('top_hits/metadata/top_hits.csv', index=False)