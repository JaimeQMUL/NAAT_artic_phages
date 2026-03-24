# This file will fetch Uniprot data, extract sequence, and extract feature table.
import requests
from typing import List, Dict
import json
from pathlib import Path

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{}.json"


def fetch_uniprot_record(accession: str) -> dict:
    """
    Fetch a UniProt record using the REST API.
    """
    response = requests.get(UNIPROT_API.format(accession))

    if response.status_code != 200:
        raise ValueError(f"UniProt accession not found: {accession}")

    return response.json()


def extract_protein_sequence(record: dict) -> str:
    """
    Extract amino acid sequence from UniProt JSON.
    """
    return record["sequence"]["value"]


def extract_feature_table(record: dict) -> List[Dict]:
    """
    Extract UniProt feature table (domains, regions, sites).
    """
    features = []

    for feature in record.get("features", []):
        if "location" in feature:
            features.append({
                "type": feature.get("type"),
                "description": feature.get("description"),
                "start": feature["location"]["start"]["value"],
                "end": feature["location"]["end"]["value"],
            })

    return features






def stage_uniprot_protein(uniprot_accession, save_directory):
    """
    Fetch UniProt protein data and store raw sequence and feature table.
    """


    # Fetch UniProt data
    record = fetch_uniprot_record(uniprot_accession)

    sequence = extract_protein_sequence(record)
    features = extract_feature_table(record)

    save_directory = Path(save_directory)

    # Save sequence
    sequence_path = save_directory / f"{uniprot_accession}_sequence.fasta"
    with open(sequence_path, "w") as f:
        f.write(f">{uniprot_accession}\n")
        f.write(sequence)

    # Save feature table
    feature_path = save_directory / f"{uniprot_accession}_features.json"
    with open(feature_path, "w") as f:
        json.dump(features, f, indent=2)

    return {
        "uniprot_accession": uniprot_accession,
        "sequence_file": str(sequence_path),
        "feature_table_file": str(feature_path),
        "num_features": len(features),
    }

stage_uniprot_protein('P04529', '../../data/uvsx')