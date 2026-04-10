import csv, json
import pandas as pd
import json
from src.biotools.fasta_tools import ExtractSequence
from pathlib import Path
from time import sleep
import requests
import xml.etree.ElementTree as ET
import numpy as np
import re


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




def WriteCleanedFasta(non_dupe_accessions, filename, protein_name):

    with open(filename, 'w') as f:
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




def BuildTaxLookup(file_path):
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    # --- CSV (NCBI style) ---
    if file_path.suffix == ".csv":
        df = pd.read_csv(file_path, low_memory=False)

        return dict(zip(df['accession'].astype(str), df['taxonId']))

    # --- JSON (UniProt style) ---
    elif file_path.suffix == ".json":
        with open(file_path) as f:
            data = json.load(f)

        return {
            entry.get('primaryAccession'): entry.get('organism', {}).get('taxonId')
            for entry in data
            if entry.get('primaryAccession')
        }

    else:
        raise ValueError(f"Unsupported file type: {file_path.suffix}")



def GetTaxonIds(df, lookups):
    taxonids = []

    for acc in df['accession']:
        acc = str(acc)
        taxid = None

        # iterate through lookup dictionaries in order
        for lookup in lookups:
            if acc in lookup:
                taxid = lookup[acc]
                break  # stop at first match

        taxonids.append(taxid)

    print(f"Taxon ids: {len(taxonids)}\nNumber of accessions: {len(df['accession'])}")

    return taxonids




def fetch_taxonomy_ranked_lineage_list(taxids, batch_size=200):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    API_KEY = "eb2c0b97f55a588259931f07f1099b896207"

    input_taxids = [str(t) for t in taxids if t is not None]

    results = []

    for i in range(0, len(input_taxids), batch_size):
        batch = input_taxids[i:i + batch_size]

        params = {
            "db": "taxonomy",
            "id": ",".join(batch),
            "retmode": "xml",
            "api_key": API_KEY
        }

        res = requests.get(url, params=params)
        root = ET.fromstring(res.text)

        # Build temporary lookup for this batch
        batch_lookup = {}

        for taxon in root.findall(".//Taxon"):
            taxid_node = taxon.find("TaxId")
            if taxid_node is None:
                continue

            taxid = taxid_node.text
            lineage_dict = {}

            lineage_ex = taxon.find("LineageEx")

            if lineage_ex is not None:
                for node in lineage_ex.findall("Taxon"):
                    rank_node = node.find("Rank")
                    name_node = node.find("ScientificName")

                    if rank_node is not None and name_node is not None:
                        rank = rank_node.text.lower()
                        name = name_node.text

                        if rank != "no rank":
                            lineage_dict[rank] = name

            current_rank = taxon.find("Rank")
            current_name = taxon.find("ScientificName")

            if current_rank is not None and current_name is not None:
                rank = current_rank.text.lower()
                name = current_name.text

                if rank != "no rank":
                    lineage_dict[rank] = name

            batch_lookup[taxid] = lineage_dict

        # ✅ Preserve original order + duplicates
        for t in batch:
            results.append({
                "taxid": t,
                "lineage": batch_lookup.get(t)
            })

        sleep(0.3)

    return results


def FindRanks(lineages, rank):
    names = []
    for i in range(len(lineages)):
        if lineages[i] is not None:
            lineage = lineages[i].get('lineage')
            if lineage is not None:
                classification = lineage.get(rank)
            else:
                classification = None
        else:
            classification = None

        names.append(classification)
    return names


def FindUniqueCounts(ranks):
    counts = {}

    for rank in ranks:
        if rank is not None:  # optional: skip missing values
            counts[rank] = counts.get(rank, 0) + 1

    return counts



def FindAccessionsInHMMResults(results_file):
    found=[]
    with open(results_file, 'r') as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith("#"):
                continue
            fields = line.split()
            accession = fields[0]
            if '.' in accession:
                accession = accession.split('.')[0]
            if '|' in accession:
                accession = accession.split('|')[1]
            found.append(accession)

        return found



def FilterFasta(fasta_file, metadata_file, output_file):


    df=pd.read_csv(metadata_file, low_memory=True)
    lengths=[]
    filtered=[]
    for acc in df['accession']:
        seq=ExtractSequence(acc,fasta_file )
        length=len(seq)
        lengths.append(length)

    mean_length = np.mean(lengths)
    std_length = np.std(lengths)

    lower = mean_length - 50
    upper = mean_length + 50
    print(upper)
    print(lower)

    for acc in df['accession']:
        seq=ExtractSequence(acc, fasta_file )
        length=len(seq)
        if length > lower and length < upper:
            filtered.append(acc)


    print(f'Upper Threshold: {upper}')
    print(f'Lower Threshold: {lower}')

    WriteCleanedFasta(filtered, output_file, 'uvsx')



def CleanNewickTree(tree):
    with open(tree) as f:
        text = f.read()

    # remove internal node labels like )0.000:
    cleaned = re.sub(r"\)\d+\.?\d*:", "):", text)

    parts = tree.split('/')[:-1]
    new_name = '/'.join(parts) + '/tree_clean.nwk'




    with open(new_name, "w") as f:
        f.write(cleaned)


def MapAccessionsToHMMResults(files):
    mapping = {}
    for file in files:
        names = file.split('/')
        name = names[-1].split('.')[0]
        found = FindAccessionsInHMMResults(file)
        print(f'{name}: Number of hits {len(found)}')
        mapping[name] = found

    return mapping




