# Script for obtaining the cifs and converting them to pdb for Molecular Dynamics

#imports
import requests, os, subprocess

bin='/opt/anaconda3/envs/cold/bin/'

# Pulling cifs from RCSB
pdb_ids= ['3IO5', '9GBG', '7Z3M']
save_dir='uvsx/data/structures/gold_standards/crystal'

for id in pdb_ids:
    url = f"https://files.rcsb.org/download/{id}.cif"

    r = requests.get(url)

    with open(f"{save_dir}/{id}.cif", "wb") as f:
        f.write(r.content)


# Convert cif to pdb
# cd uvsx/data/structures/gold_standards/crystal
# maxit -input 3IO5.cif -output 3IO5.pdb -o 2
# maxit -input 9GBG.cif -output 9GBG.pdb -o 2
# maxit -input 7Z3M.cif -output 7Z3M.pdb -o 2






