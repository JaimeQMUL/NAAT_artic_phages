import requests
import json

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
    with open(f"{protein_name}data/curated_database/interpro_domain_matches_metadata.json") as f:
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


    #I find 1048 instead of 1302
