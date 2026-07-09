import os, subprocess
from src.biotools.fasta_tools import *

# Create HMM profiles and use them to search database

# Add new line at the end of every fasta file
os.chdir('../1-DatabaseCuration/gold_standards/')
subprocess.run('cd references && for f in *.fasta; do tail -c 1 "$f" | read -r _ || printf "\n" >> "$f"; done', shell=True)

# =================================================== GOLD STANDARDS HMM ===============================================
# First HMM is of the full alignment of validated seqs. ['7Z3M.fasta', '9GBG.fasta', 'K35G_E198R.fasta', 'P04529.fasta']
# Combining gold standards into one fasta
subprocess.run('cat P04529.fasta K35G_E198R.fasta 7Z3M.fasta 9GBG.fasta > gold_standards.fasta', shell=True)

#Aligning gold standards using famsa
subprocess.run('/opt/anaconda3/envs/cold/bin/famsa gold_standards.fasta ../2-HiddenMarkovModels/seeds/gold_standards.aln.fasta', shell=True)

#Convert fasta to stockholm
os.chdir('../2-HiddenMarkovModels/seeds/')
subprocess.run('/opt/anaconda3/envs/cold/bin/esl-reformat stockholm gold_standards.aln.fasta > gold_standards.aln.sto', shell=True)
subprocess.run('rm gold_standards.aln.fasta', shell=True)

os.chdir('..')
subprocess.run('/opt/anaconda3/envs/cold/bin/hmmbuild profiles/gold_standards.hmm seeds/gold_standards.aln.sto', shell=True)

# ==================================================== PFAM HMMs =======================================================
# HMMS obtained from interpro using ids
interpro_ids=['PF21134', 'PF00154']
for id in interpro_ids:
    # download
    subprocess.run(
        f'curl -L "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{id}/?annotation=alignment:seed&download" '
        f'--compressed -o seeds/{id}_seed.sto.gz',
        shell=True,
        check=True
    )

    # unzip
    subprocess.run(
        f'gunzip -f seeds/{id}_seed.sto.gz',
        shell=True,
        check=True
    )

    # Build profile from seeds
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/hmmbuild profiles/{id}.hmm seeds/{id}_seed.sto', shell=True)


# ================================================= CUSTOM DOMAIN HMM ==================================================
# Will use interpro id HMMs to scan reference fasta for domains, extract coords of domains found, align them and use them
# as a seed for hmm profiles.

for id in interpro_ids:
    #prepare hmm database
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/hmmpress profiles/{id}.hmm', shell=True)
    # Search gold standards for coords of domain matches
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/hmmscan --domtblout tmp/{id}_results.domtblout profiles/{id}.hmm ../1-DatabaseCuration/gold_standards/gold_standards.fasta', shell=True)
    #remove index files
    subprocess.run('rm profiles/*.hmm.h*', shell=True)

    #extract coords from hmmscan and index protein and store domain sequence in a fasta
    CreateDomainFasta(f'tmp/{id}_results.domtblout', f'tmp/{id}_domain.fasta', '../1-DatabaseCuration/gold_standards/P04529.fasta')

    # Aligning domains
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/famsa tmp/{id}_domain.fasta seeds/{id}_domain.aln.fasta', shell=True)

    # Converting to stockholm
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/esl-reformat stockholm seeds/{id}_domain.aln.fasta > seeds/{id}_domain.aln.sto', shell=True)
    subprocess.run(f'rm seeds/{id}_domain.aln.fasta', shell=True)

    # Create hmm profile from alignment
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/hmmbuild profiles/{id}_domain.hmm seeds/{id}_domain.aln.sto', shell=True)


# ====================================== Searching curated database using all 5 HMMs ===================================
profiles=os.listdir('profiles/')
for profile in profiles:
    name=profile.split('.')[0]
    subprocess.run(f'/opt/anaconda3/envs/cold/bin/hmmsearch --tblout results/{name}_results.tblout profiles/{profile} ../1-DatabaseCuration/clean/sequences/cleaned_curated_database.fasta', shell=True)

# Summary of results
for result in os.listdir('results/'):
    count=0
    with open(f'results/{result}') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            count+=1
        print(f'{result}\t Number of hits: {count}')


# cleanup
# remove files in tmp
subprocess.run('rm tmp/*', shell=True)