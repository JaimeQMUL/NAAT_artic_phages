

import requests
import json
from time import sleep
import csv



########################################################################################################################
# Sourcing all proteins from uniprot using search queries
########################################################################################################################
def QuerySearchUniprot(query, protein_name):

    # UniProt REST API endpoint for JSON
    url = "https://rest.uniprot.org/uniprotkb/search"

    # Start params
    params = {
        "query": query,
        "format": "json",  # JSON returns *all available metadata*
        "size": 500  # max per request; use pagination for more
    }

    all_results = []

    while True:
        response = requests.get(url, params=params)
        if response.status_code != 200:
            raise Exception(f"API request failed with status code {response.status_code}")

        data = response.json()
        all_results.extend(data.get("results", []))  # append entries

        # Check for pagination
        next_link = None
        if "Link" in response.headers:
            links = response.headers["Link"].split(",")
            for l in links:
                if 'rel="next"' in l:
                    next_link = l.split(";")[0].strip("<>")
                    break
        if next_link:
            url = next_link
            params = {}  # already included in next_link
        else:
            break

    # Save all metadata as JSON
    with open(f"{protein_name}/data/curated_database/uniprot_query_search_metadata.json", "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"Downloaded {len(all_results)} entries with full metadata")


    # Save sequences extracted in metadata to fasta file

    # Load your JSON metadata file from UniProt
    with open(f"{protein_name}/data/curated_database/uniprot_query_search_metadata.json") as f:
        data = json.load(f)

    # Output FASTA file
    with open(f"{protein_name}/data/curated_database/uniprot_query_search.fasta", "w") as fasta_file:
        for entry in data:
            accession = entry.get("primaryAccession")  # UniProt accession
            protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")
            sequence = entry.get("sequence", {}).get("value")

            if sequence:  # make sure there is a sequence
                header = f">{accession} {protein_name} [{organism}]"
                fasta_file.write(header + "\n")

                # Wrap sequence at 60 characters per line (optional, standard FASTA formatting)
                for i in range(0, len(sequence), 60):
                    fasta_file.write(sequence[i:i+60] + "\n")

    print(f"Saved {len(data)} uniprot_query_search.fasta")


########################################################################################################################
# Sourcing all proteins from ncbi using search queries
########################################################################################################################

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


########################################################################################################################
# Sourcing all proteins that contain a InterPro domain id match from Uniprot
########################################################################################################################


def InterproSearchUniprot(interpro_ids, protein_name):
    # Build the query
    query = " OR ".join([f"(database:InterPro {ipr})" for ipr in interpro_ids])

    # UniProt API URL (TSV format with only accession)
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "json",
        "size": 500             # max per page
    }

    all_results = []
    next_url = f"{base_url}"

    while next_url:
        r = requests.get(next_url, params=params)
        if r.status_code != 200:
            print("Error:", r.status_code)
            break

        data = r.json()

        for entry in data["results"]:
            lineage = entry.get('organism').get('lineage')
            if lineage[0] == 'Viruses':
                all_results.append(entry)

        # Check for next page
        next_url = None
        if "Link" in r.headers:
            for link in r.headers["Link"].split(","):
                if 'rel="next"' in link:
                    next_url = link.split(";")[0].strip("<>")

        # After first request, params are in the URL for pagination
        params = None

    # Save all metadata as JSON
    with open(f"{protein_name}/data/curated_database/interpro_domain_matches_metadata.json", "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"Downloaded {len(all_results)} entries with full metadata")

    # Load your JSON metadata file from UniProt
    with open(f"{protein_name}/data/curated_database/interpro_domain_matches_metadata.json") as f:
        data = json.load(f)

    # Output FASTA file
    with open(f"{protein_name}/data/curated_database/interpro_domain_matches.fasta", "w") as fasta_file:
        for entry in data:
            accession = entry.get("primaryAccession")  # UniProt accession
            protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")
            sequence = entry.get("sequence", {}).get("value")

            if sequence:  # make sure there is a sequence
                header = f">{accession} {protein_name} [{organism}]"
                fasta_file.write(header + "\n")

                # Wrap sequence at 60 characters per line (optional, standard FASTA formatting)
                for i in range(0, len(sequence), 60):
                    fasta_file.write(sequence[i:i+60] + "\n")

    print(f"Saved {len(data)} sequences to interpro_domain_matches.fasta")



########################################################################################################################
# Sourcing proteins stored in the RCSB database using IDs found in the literature
########################################################################################################################

BASE = "https://data.rcsb.org/rest/v1/core"


def get_entry(pdb_id):
    return requests.get(f"{BASE}/entry/{pdb_id}").json()


def get_polymer_entity(pdb_id, entity_id):
    return requests.get(f"{BASE}/polymer_entity/{pdb_id}/{entity_id}").json()


def IDSearchRCSB(pdb_ids, protein_name):
    for pdb_id in pdb_ids:
        entry = get_entry(pdb_id)

        # --- Save metadata (everything except sequence) ---
        with open(f"{protein_name}/data/references/{pdb_id}_metadata.json", "w") as f:
            json.dump(entry, f, indent=2)

        # --- Extract sequences ---
        entity_ids = entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]

        fasta_lines = []

        for entity_id in entity_ids:
            entity = get_polymer_entity(pdb_id, entity_id)

            sequence = entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
            description = entity["rcsb_polymer_entity"]["pdbx_description"]

            header = f">{pdb_id}_entity_{entity_id} {description}"
            fasta_lines.append(header)
            fasta_lines.append(sequence.replace("\n", ""))  # clean formatting

        # --- Save FASTA ---
        with open(f"{protein_name}/data/references/{pdb_id}.fasta", "w") as f:
            f.write("\n".join(fasta_lines))




########################################################################################################################
# Recreates identified mutations that improve performance identified through wet lab research
########################################################################################################################

# Single Amino acid mutations

def CreateMutants(mutants, reference, protein_name):



    positions = {}

    for mutant in mutants:
        parts = mutant.split('/')  # handles both single + double

        pos_list = []
        for part in parts:
            pos = int(part[1:-1]) - 1
            pos_list.append(pos)

        positions[mutant] = pos_list

    for item in positions.items():
        new_base=item[0][-1]
        name=item[0]
        reference_list=list(reference)
        for pos in item[1]:
            reference_list[pos]=new_base
        new_seq=''.join(reference_list)
        safe_name = name.replace("/", "_")
        with open(f"{protein_name}/data/references/{safe_name}.fasta",'w') as f:
            f.write(f">{safe_name} Rational Design UvsX mutant\n{new_seq}")



##########################
# Larger scale mutations #
##########################
# # Getting sequence of loop 2 donor recombinase
# # IN Future develop this to recreate published mutations involving larger scale modifications, in this case a full loop swap.
# header_to_find = ">YP_003097304"
# sequence_lines = []
# found = False
#
# with open('../data/ncbi/uvsx_sequences.fasta', 'r') as f:
#     for line in f:
#         line = line.strip()
#         if line.startswith(">"):
#             if found:
#                 # We've reached the next header, stop reading
#                 break
#             if line.startswith(header_to_find):
#                 found = True
#         elif found:
#             # Collect sequence lines
#             sequence_lines.append(line)
#
# # Combine all sequence lines into a single string
# # full_sequence = "".join(sequence_lines)
# # print(f"{header_to_find}:\n{full_sequence}")
# # idea for this was replicating the loop 2 swap from YP_003097304 into UvsX
