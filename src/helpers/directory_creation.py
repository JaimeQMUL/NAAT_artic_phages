import os

def CreateDirectoryStructure(protein_name):
    os.makedirs(f"{protein_name}/data", exist_ok=True)
    os.makedirs(f"{protein_name}/data/curated_database", exist_ok=True)
    os.makedirs(f"{protein_name}/data/references", exist_ok=True)