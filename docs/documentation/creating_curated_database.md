# UvsX Database Generation

## Summary

This script builds a dataset of UvsX-like protein sequences by querying multiple biological databases, extracting sequences and metadata, and generating additional reference variants.

---

## Steps Performed

### 1. Reference Sequence Retrieval
- Downloaded the UvsX reference protein sequence from UniProt  
- Accession: `P04529`  
- Saved as a FASTA file in the `references` directory  

---

### 2. Directory Setup
- Created a structured directory system to store:
  - Raw sequences
  - Metadata
  - Reference data

---

### 3. UniProt Search
- Queried UniProt using:
(recA OR uvsX OR recombinase) AND taxonomy_id:10239
- Retrieved:
- Full metadata (JSON)
- Protein sequences (FASTA)
- Used pagination to collect all results

---

### 4. NCBI Search
- Queried NCBI Protein database using:
(uvsx OR recA) AND txid10239[Organism:exp]

- Retrieved:
- Protein sequences (FASTA)
- Metadata (CSV)
- Used batching to handle large result sets

---

### 5. InterPro Domain Search
- Queried UniProt for proteins containing InterPro domains:
- `IPR049428`
- `IPR049047`
- Filtered results to viral proteins only
- Saved:
- Metadata (JSON)
- Sequences (FASTA)

---

### 6. Mutant Sequence Generation
- Generated mutant sequences from literature-reported mutations:
- `E198N`, `E198R`, `E198K`, `K35G`, `K35G/E198R`, `D274A`
- Applied mutations to the reference sequence
- Saved each mutant as a FASTA file

---

### 7. RCSB PDB Data Extraction
- Retrieved protein data for selected PDB IDs:
- `7Z3M`
- `9GBG`
- Saved:
- Full metadata (JSON)
- Extracted sequences (FASTA)

---

## Output

The pipeline produces:

- FASTA files containing protein sequences from:
- UniProt
- NCBI
- InterPro searches
- Mutants
- PDB entries

- Metadata files:
- JSON (UniProt, InterPro, PDB)
- CSV (NCBI)

---

## Notes

- API pagination and batching were used to retrieve large datasets
- Delays (`sleep`) were added to avoid rate limiting
- All sequences were standardised into FASTA format for downstream analysis
