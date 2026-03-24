import requests
import json
########################################################################################################################
# Sourcing all UvsX/RecA-like proteins from search queries
########################################################################################################################

# Broad query for viral RecA/UvsX/recombinase proteins
query = "(recA OR uvsX OR recombinase) AND taxonomy_id:10239"

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
with open("../../data/uniprot/viral_recA_uvsX_all_metadata.json", "w") as f:
    json.dump(all_results, f, indent=2)

print(f"Downloaded {len(all_results)} entries with full metadata")


# Save sequences extracted in metadata to fasta file

# Load your JSON metadata file from UniProt
with open("../../data/uniprot/viral_recA_uvsX_all_metadata.json") as f:
    data = json.load(f)

# Output FASTA file
with open("../../data/uniprot/viral_recA_uvsX.fasta", "w") as fasta_file:
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

print(f"Saved {len(data)} sequences to viral_recA_uvsX.fasta")