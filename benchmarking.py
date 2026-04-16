import subprocess
import pandas as pd
from src.biotools.fasta_tools import *
from time import time
from src.qc.fasta_qc import *
import os
# BENCHMARKING DIFFERENT TOOLS I AM USING

os.chdir("uvsx/data/tests/")

# ALIGNMENTS
##################################### Creating subset of the full filtered fasta #######################################

df=pd.read_csv('../curated_database/filtered_metadata.csv')

df_sample = df.sample(n=100, random_state=0)

# Add UvsX t4 to data
uvsx_t4=df[df['accession']=='P04529']
df_sample = pd.concat([df_sample, uvsx_t4], ignore_index=True)

total_lengths = 0

total_lengths=0


with open('subset_filtered_database.fasta', 'w') as f:
    for acc in df_sample['accession']:
        header, seq=ExtractSequence(acc, '../curated_database/filtered_curated_database.fasta', header=True)
        f.write(f'{header}\n{seq}\n')
        total_lengths+=len(seq)

# calculating mean length
mean_length=total_lengths/len(df_sample)
########################################################################################################################

######################################################## MAFFT #########################################################
# Default parameters

# MAFFT alignment
start = time()
subprocess.run('/opt/anaconda3/envs/cold/bin/mafft subset_filtered_database.fasta > subset_filtered_database_default_alignment.fasta', shell=True)
end = time()

#Trim gaps using Trimal
subprocess.run('/opt/anaconda3/envs/cold/bin/trimal -automated1 -in subset_filtered_database_default_alignment.fasta -out subset_filtered_database_default_alignment_trimmed.fasta', shell=True)

print('#'*50, 'Default MAFFT', '#'*50)
print(f'Total time: {end-start}')
print(f'Default MAFFT alignment took {round(end-start,4)} seconds')
mean_gap, gap_counts = CountAlignmentGaps('subset_filtered_database_default_alignment.fasta')
trim_mean_gap, trim_gap_counts = CountAlignmentGaps('subset_filtered_database_default_alignment_trimmed.fasta')
print(f'Average sequence Length is {mean_length}\nAfter Alignment Avergage Gap count per sequence is {mean_gap}\nAfter Trimming Average Gap count per sequence is {trim_mean_gap}')

# Auto parameters
start = time()
# MAFFT alignment
subprocess.run('/opt/anaconda3/envs/cold/bin/mafft --auto subset_filtered_database.fasta > subset_filtered_database_auto_alignment.fasta', shell=True)
end = time()
#Trim gaps using Trimal
subprocess.run('/opt/anaconda3/envs/cold/bin/trimal -automated1 -in subset_filtered_database_auto_alignment.fasta -out subset_filtered_database_auto_alignment_trimmed.fasta', shell=True)

print('#'*50, 'Auto MAFFT', '#'*50)
print(f'Total time: {end-start}')
print(f'Auto MAFFT alignment took {round(end-start,4)} seconds')
mean_gap, gap_counts = CountAlignmentGaps('subset_filtered_database_auto_alignment.fasta')
trim_mean_gap, trim_gap_counts = CountAlignmentGaps('subset_filtered_database_auto_alignment_trimmed.fasta')
print(f'Average sequence Length is {mean_length}\nAfter Alignment Avergage Gap count per sequence is {mean_gap}\nAfter Trimming Average Gap count per sequence is {trim_mean_gap}')


# L-Ins-I alignment
start = time()
subprocess.run('/opt/anaconda3/envs/cold/bin/linsi subset_filtered_database.fasta > subset_filtered_database_linsi_alignment.fasta', shell=True)
end = time()

#Trim gaps using Trimal
subprocess.run('/opt/anaconda3/envs/cold/bin/trimal -automated1 -in subset_filtered_database_linsi_alignment.fasta -out subset_filtered_database_linsi_alignment_trimmed.fasta',shell=True)

print('#'*50, 'Linsi MAFFT', '#'*50)
print(f'Total time: {end-start}')
print(f'Linsi MAFFT alignment took {round(end-start,4)} seconds')
mean_gap, gap_counts = CountAlignmentGaps('subset_filtered_database_linsi_alignment.fasta')
trim_mean_gap, trim_gap_counts = CountAlignmentGaps('subset_filtered_database_linsi_alignment.fasta')
print(f'Average sequence Length is {mean_length}\nAfter Alignment Average Gap count per sequence is {mean_gap}\nAfter Trimming Average Gap count per sequence is {trim_mean_gap}')

######################################################## FAMSA #########################################################
# Alignment
start = time()
subprocess.run('/opt/anaconda3/envs/cold/bin/famsa subset_filtered_database.fasta subset_filtered_database_famsa_alignment.fasta', shell=True)
end = time()

# Alignment quality
print('#'*50, 'Default FAMSA', '#'*50)
print(f'Total time: {end-start}')
mean_gap, gap_counts = CountAlignmentGaps('subset_filtered_database_famsa_alignment.fasta')
mean_occ, occ_counts = AlignmentOccupancy('subset_filtered_database_famsa_alignment.fasta')

#Trim alignment
subprocess.run('/opt/anaconda3/envs/cold/bin/trimal -automated1 -in subset_filtered_database_famsa_alignment.fasta -out subset_filtered_database_famsa_alignment_trimmed.fasta', shell=True)
trim_mean_gap, trim_gap_counts = CountAlignmentGaps('subset_filtered_database_famsa_alignment_trimmed.fasta')
trim_mean_occ, trim_occ_counts = AlignmentOccupancy('subset_filtered_database_famsa_alignment_trimmed.fasta')
print(f'Average sequence Length is {mean_length}\nAfter Alignment Average Gap count per sequence is {mean_gap}\nAfter Trimming Average Gap count per sequence is {trim_mean_gap}')
print(f'Mean Occupancy before Trim: {mean_occ}\nMean Occupancy after Trim {trim_mean_occ}')


# Will proceed with famsa for alignment as it is much faster for


# alignment tree
# FAMSA TREE
subprocess.run('/opt/anaconda3/envs/cold/bin/famsa -gt nj -gt_export subset_filtered_database.fasta famsa_tree.dnd', shell=True)

# IQtree
subprocess.run("/opt/anaconda3/envs/cold/bin/iqtree3 -s subset_filtered_database_famsa_alignment_trimmed.fasta", shell=True)