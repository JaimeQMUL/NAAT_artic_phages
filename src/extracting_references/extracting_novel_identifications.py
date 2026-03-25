import requests
import json

# Sequences identified in paper from ice samples. 7Z3M is thermophilic and 9GBG is psychrophilic  https://academic.oup.com/nar/article/54/4/gkag069/8462617?login=true

pdb_ids=['7Z3M','9GBG']

BASE = "https://data.rcsb.org/rest/v1/core"


def get_entry(pdb_id):
    return requests.get(f"{BASE}/entry/{pdb_id}").json()


def get_polymer_entity(pdb_id, entity_id):
    return requests.get(f"{BASE}/polymer_entity/{pdb_id}/{entity_id}").json()


for pdb_id in pdb_ids:
    entry = get_entry(pdb_id)

    # --- Save metadata (everything except sequence) ---
    with open(f"../../data/uvsx/{pdb_id}_metadata.json", "w") as f:
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
    with open(f"../../data/uvsx/{pdb_id}.fasta", "w") as f:
        f.write("\n".join(fasta_lines))