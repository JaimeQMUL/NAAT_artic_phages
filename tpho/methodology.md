# 🧪 Methodology: Phosphatase Sequence Analysis and Model Preparation

## 1. Data Preparation

1. Collect:
   - The **predicted phosphatase sequence** (target protein)
   - The **training dataset sequences**

2. Combine all sequences into a single FASTA file:
   - Create a dedicated folder (e.g., `FASTA/`)
   - Save all sequences in standard FASTA format within this folder

---

## 2. Multiple Sequence Alignment and Phylogenetic Analysis

1. Perform multiple sequence alignment using ClustalW:
   - https://www.genome.jp/tools-bin/clustalw  

2. Download the resulting:
   - Alignment output
   - **DND (guide tree) file**

3. Import the DND file into MEGA11

4. Construct and visualise the **phylogenetic tree**

5. Identify:
   - The position of the **target phosphatase (UniProt ID)** within the tree  
   - Closely related sequences

6. Based on clustering, assign or infer the **enzyme EC number** for the target phosphatase

---

## 3. Dataset Annotation

1. Update the dataset file:
   - `original_Complete_Information_In_sequence.tsv`

2. Add:
   - The inferred **EC number**
   - Relevant phosphatase annotation information

---

## 4. EC-Based Sequence Grouping

1. Extract:
   - The target phosphatase sequence  
   - All sequences sharing the same EC number  

2. Save these sequences into a FASTA file:
   - Store in the `ComparePhosphatase/` directory

---

## 5. Multiple Sequence Alignment (Local Processing)

1. Run the script:
   - `MultipleSequenceAlignment.py`

2. Input:
   - FASTA file from Step 4

3. Output:
   - `results.tsv` containing:
     - Conserved amino acid positions
     - Alignment-derived conservation indices

---

## 6. Intermediate Data Integration

1. Transfer relevant alignment results into:
   - `sequence_Alignment_Results_Intermediate_File.tsv`

---

## 7. Sequence Processing (Removal of Conserved Residues)

1. Run:
   - `ProcessingComparedSequence.py`

2. Function:
   - Removes **conserved amino acids** from sequences

3. Output:
   - `test_data.tsv` (processed sequence data)

---

## 8. FASTA Generation for Model Input

1. Convert processed sequences into FASTA format:
   - File: `test_data.fasta`

2. Ensure:
   - Conserved residues are removed
   - Sequence IDs remain consistent

---

## 9. Model Execution

1. Run:
   - `test_model.py`

2. Input:
   - `test_data.fasta`

3. Purpose:
   - Perform prediction or classification on processed phosphatase sequences

---

## Workflow Summary
FASTA prep
↓
ClustalW alignment → DND tree
↓
MEGA11 phylogeny → EC assignment
↓
Dataset annotation (TSV)
↓
EC-based grouping
↓
Local MSA script
↓
Extract conserved residues
↓
Remove conserved residues
↓
Generate test FASTA
↓
Run model