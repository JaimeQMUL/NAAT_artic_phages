import csv, json


def GetUniqueAccessions(protein_name):

    #Find the accession from metadata
    ncbi_accessions = []
    csv_path = f"../{protein_name}/data/curated_database/ncbi_query_search_metadata.csv"
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            accession = row.get("accession")
            if accession:
                ncbi_accessions.append(accession)
    print(f"NCBI Accessions: {len(ncbi_accessions)}")

    # Get all query search accessions
    uniprot_search_accessions = []
    with open(f"../{protein_name}/data/curated_database/uniprot_query_search_metadata.json") as f:
        data = json.load(f)
        for entry in data:
            accession = entry.get("primaryAccession")  # UniProt accession
            uniprot_search_accessions.append(accession)
    print(f"Uniprot Accessions: {len(uniprot_search_accessions)}")

    # Get all interpro match accessions
    domain_match_accessions = []
    with open(f"../{protein_name}/data/curated_database/interpro_domain_matches_metadata.json") as f:
        data = json.load(f)
        for entry in data:
            accession = entry.get("primaryAccession")  # UniProt accession
            domain_match_accessions.append(accession)
    print(f"Domain match Accessions: {len(domain_match_accessions)}")

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

