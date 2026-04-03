import csv, json
import pandas as pd
import json
from src.biotools.fasta_tools import ExtractSequence


def GetUniqueAccessions(protein_name):

    #Find the accession from metadata
    ncbi_accessions = []
    csv_path = f"{protein_name}/data/curated_database/ncbi_query_search_metadata.csv"
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            accession = row.get("accession")
            if accession:
                ncbi_accessions.append(accession)
    print(f"NCBI Accessions: {len(ncbi_accessions)}")

    # Get all query search accessions
    uniprot_search_accessions = []
    with open(f"{protein_name}/data/curated_database/uniprot_query_search_metadata.json") as f:
        data = json.load(f)
        for entry in data:
            accession = entry.get("primaryAccession")  # UniProt accession
            uniprot_search_accessions.append(accession)
    print(f"Uniprot Accessions: {len(uniprot_search_accessions)}")

    # Get all interpro match accessions
    domain_match_accessions = []
    with open(f"{protein_name}/data/curated_database/interpro_domain_matches_metadata.json") as f:
        data = json.load(f)
        for entry in data:
            accession = entry.get("primaryAccession")  # UniProt accession
            domain_match_accessions.append(accession)
    print(f"Domain match Accessions: {len(domain_match_accessions)}")

    # Create lists of accession from each file
    unique_accessions = []
    for acc in ncbi_accessions:
        if acc not in unique_accessions:
            unique_accessions.append(acc)
    for acc in uniprot_search_accessions:
        if acc not in unique_accessions:
            unique_accessions.append(acc)
    for acc in domain_match_accessions:
        if acc not in unique_accessions:
            unique_accessions.append(acc)

    return unique_accessions


# Full script

from collections import defaultdict


def RemoveRedundancies(unique_accessions, protein_name):
    files = [
        f'{protein_name}/data/curated_database/uniprot_query_search.fasta',
        f'{protein_name}/data/curated_database/ncbi_query_search.fasta',
        f'{protein_name}/data/curated_database/interpro_domain_matches.fasta'
    ]

    # NCBI CSV
    ncbi_df = pd.read_csv(f'{protein_name}/data/curated_database/ncbi_query_search_metadata.csv', low_memory=False)

    # UniProt JSON
    with open(f'{protein_name}/data/curated_database/uniprot_query_search_metadata.json') as f:
        uniprot_data = json.load(f)

    # InterPro JSON
    with open(f'{protein_name}/data/curated_database/interpro_domain_matches_metadata.json') as f:
        interpro_data = json.load(f)

    ncbi_lookup = dict(zip(ncbi_df["accession"], ncbi_df["taxonId"]))

    uniprot_lookup = {
        entry.get("primaryAccession"): entry.get("organism", {}).get("taxonId")
        for entry in uniprot_data
    }

    interpro_lookup = {
        entry.get("primaryAccession"): entry.get("organism", {}).get("taxonId")
        for entry in interpro_data
    }

    seq_map = {}

    #Map accessions to sequences, will show overlaps in accessions
    for accession in unique_accessions:
        seq = ''
        for file in files:
            seq = ExtractSequence(accession, file)
            if seq:
                break
        seq_map[accession] = seq

    seq_tax_groups = defaultdict(list)

    for acc, seq in seq_map.items():
        if not seq:
            continue

        tax = (
                uniprot_lookup.get(acc)
                or ncbi_lookup.get(acc)
                or interpro_lookup.get(acc)
        )

        tax = tax or "UNKNOWN"

        key = (seq, tax)  # <-- critical change
        seq_tax_groups[key].append(acc)

    duplicates = {
        key: accs for key, accs in seq_tax_groups.items() if len(accs) > 1
    }

    non_dupe_accessions = []

    for (seq, tax), accs in seq_tax_groups.items():
        # keep one representative accession per (sequence, taxon)
        non_dupe_accessions.append(accs[0])


    print(len(non_dupe_accessions))
    return non_dupe_accessions




def WriteCleanedFasta(non_dupe_accessions, protein_name):

    with open(f'{protein_name}/data/curated_database/cleaned_curated_database.fasta', 'w') as f:
        for acc in non_dupe_accessions:
            header, seq = ExtractSequence(acc, f'{protein_name}/data/curated_database/ncbi_query_search.fasta', header=True)
            if seq=='':
                header, seq = ExtractSequence(acc, f'{protein_name}/data/curated_database/uniprot_query_search.fasta', header=True)
            f.write(f'{header} \n')
            for i in range(0, len(seq), 60):
                f.write(f'{seq[i:i + 60]}\n')




def FindAccession(accession, file_path):
    accession = str(accession)

    # -------------------------
    # CSV (e.g. NCBI metadata)
    # -------------------------
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path, low_memory=False)

        # Try common column names
        possible_cols = ['accession', 'caption', 'Accession', 'primaryAccession']

        for col in possible_cols:
            if col in df.columns:
                match = df[df[col].astype(str) == accession]
                if not match.empty:
                    return match.to_dict(orient='records')[0]

        return None

    # -------------------------
    # JSON (e.g. UniProt / NCBI)
    # -------------------------
    elif file_path.endswith(".json"):
        with open(file_path) as f:
            data = json.load(f)
            for entry in data:
                if entry.get("primaryAccession") == accession:
                    return entry


        return None

    else:
        raise ValueError("Unsupported file type (must be .csv or .json)")




