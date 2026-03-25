# Extracting the known UvsX-like sequences from NCBI as well as metadata.

import requests
from time import sleep
import csv

def QuerySearchNCBI(query, protein_name):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    all_ids = []
    retstart = 0
    retmax = 1000

    while True:
        params = {
            "db": "protein",
            "term": query,
            "retstart": retstart,
            "retmax": retmax,
            "retmode": "json"
        }

        res = requests.get(url, params=params).json()
        ids = res["esearchresult"]["idlist"]

        if not ids:
            break

        all_ids.extend(ids)
        retstart += retmax

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "protein",
        "id": ",".join(all_ids[:200]),  # batch of IDs
        "retmode": "json"
    }

    res = requests.get(url, params=params, verify=True)
    data = res.json()

    # extract accessions
    accessions = []

    for uid, rec in data["result"].items():
        if uid == "uids":
            continue  # skip the list of IDs
        accessions.append(rec.get("caption"))


    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    batch_size = 200
    all_sequences = ""

    for start in range(0, len(all_ids), batch_size):
        batch = all_ids[start:start + batch_size]
        params = {
            "db": "protein",
            "id": ",".join(batch),
            "rettype": "fasta",
            "retmode": "text"
        }
        res = requests.get(fetch_url, params=params)
        all_sequences += res.text
        print(f"Fetched sequences {start + 1} to {start + len(batch)}")
        sleep(0.4)  # polite delay to avoid hitting NCBI limits

        with open(f"{protein_name}/data/curated_database/ncbi_query_search.fasta", "w") as f:
            f.write(all_sequences)

    print(f"All sequences saved to '{protein_name}/data/curated_database/ncbi_query_search.fasta'")


    batch_size = 200
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    metadata_records = []

    for start in range(0, len(all_ids), batch_size):
        batch = all_ids[start:start + batch_size]
        params = {
            "db": "protein",
            "id": ",".join(batch),
            "retmode": "json"
        }
        res = requests.get(summary_url, params=params).json()

        # Process this batch immediately
        for uid, rec in res["result"].items():
            if uid == "uids":  # skip the uids list
                continue

            subtype = rec.get("subtype")
            subname = rec.get("subname")

            # Base record
            record = {
                "uid": rec.get("uid"),
                "accession": rec.get("caption"),
                "title": rec.get("title"),
                "organism": rec.get("organism"),
                "sequence_length": rec.get("slen"),
                "create_date": rec.get("createdate"),
                "update_date": rec.get("updatedate"),
                "dbsource": rec.get("sourcedb"),
                "extra": rec.get("extra")
            }

            # Add subtype/subname if present
            if subtype and subname:
                type_parts = subtype.split('|')
                name_parts = subname.split('|')
                for t, n in zip(type_parts, name_parts):
                    if t and n:
                        record[t] = n

            metadata_records.append(record)

        print(f"Processed batch {start+1} to {start+len(batch)}")
        sleep(0.4)

    print(f"Total records collected: {len(metadata_records)}")

    with open(f"{protein_name}/data/curated_database/ncbi_query_search_metadata.csv", "w", newline="", encoding="utf-8") as csvfile:
        fieldnames=[]
        for i in metadata_records:
            fields=i.keys()
            for f in fields:
                if f not in fieldnames:
                    fieldnames.append(f)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metadata_records)

    print(f"Metadata saved to '{protein_name}/data/curated_database/ncbi_query_search_metadata.csv'")