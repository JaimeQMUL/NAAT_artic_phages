# Results

## 1. Identification of UvsX orthologues from sequence databases

### 1.1 Construction of the curated UvsX sequence dataset
- Show the number of sequences retrieved from:
  - UniProt
  - NCBI
  - InterPro
- Show reduction in sequence numbers after:
  - duplicate accession removal
  - identical sequence/species filtering
  - sequence length filtering
- Present the final curated dataset size

**Figures/Tables:**
- Database workflow figure
- Filtering summary table

---

### 1.2 HMM-based identification of UvsX candidates
- Present performance/output of the three HMM searches
- Show:
  - number of hits from each HMM
  - overlap between HMMs
  - sequences retained after profile screening

**Figures/Tables:**
- Venn diagram/upset plot of HMM hits
- HMM hit statistics table

---

### 1.3 Motif-based validation of candidate UvsX proteins
- Present Walker A motif detection results
- Present Walker B motif detection results
- Show how motif filtering reduced the candidate pool
- Explain final five-criterion filtering approach

**Figures/Tables:**
- Filtering funnel
- Summary table of screening criteria

---

### 1.4 Final high-confidence UvsX candidate dataset
- Present the final 1,498 candidates
- Describe:
  - taxonomic distribution
  - host organism diversity
  - sequence diversity
  - similarity to gold-standard proteins

**Figures/Tables:**
- Taxonomic distribution plot
- Sequence similarity distribution
- Candidate summary table

---

# 2. Thermal prediction and classification of candidate UvsX proteins

## 2.1 Seq2Topt predicted thermal optima
- Present predicted temperature distribution across all candidates
- Identify:
  - lowest predicted Topt
  - highest predicted Topt
  - average/median Topt
- Describe whether candidates are enriched towards specific temperature ranges

**Figures/Tables:**
- Histogram/density plot of predicted Topt
- Summary statistics table

---

## 2.2 Thermal classification using TMP
- Introduce benchmark dataset
- Present composition of validation dataset:
  - psychrophilic
  - mesophilic
  - thermophilic examples

**Figures/Tables:**
- Benchmark dataset table

---

## 2.3 Comparison of TMP implementations
Compare:

### Fast TMP
- Accuracy
- Strengths/limitations

### Democratic TMP
- Accuracy
- Effect of multiple evidence sources

### Summary TMP
- Accuracy
- Improvement over previous approaches

**Figures/Tables:**
- Accuracy comparison plot
- Confusion matrices
- Performance metrics table

---

## 2.4 Agreement between TMP and Seq2Topt
- Compare literature-derived classifications with sequence predictions
- Identify:
  - agreeing candidates
  - conflicting candidates
- Discuss whether combining approaches improves confidence

**Figures/Tables:**
- Scatter plot of Seq2Topt temperature vs TMP class
- Agreement matrix

---

# 3. Evolutionary analysis of UvsX candidates (optional)

*(Only include if you use this later)*

## 3.1 Phylogenetic placement of candidate UvsX proteins
- Show relationship between:
  - gold-standard proteins
  - candidate proteins
  - thermal groups

## 3.2 Relationship between phylogeny and thermal adaptation
- Examine whether thermal groups cluster
- Identify evolutionary patterns

**Figures/Tables:**
- Phylogenetic tree
- iTOL visualisation

---

# 4. Structural prediction and validation of UvsX proteins

## 4.1 Comparison of structure prediction methods
Compare:

- Boltz
- AlphaFold

Using:
- crystal structures
- US-align scores
- RMSD/TM-score

**Figures/Tables:**
- Structure overlays
- Model comparison table

---

## 4.2 Effect of sequence input strategy on Boltz prediction
Compare:

- Full-length sequences
- Trimmed sequences
- UniRef100 MSA input
- Custom UvsX MSA input

Determine:
- best prediction strategy
- whether custom alignments improve accuracy

**Figures/Tables:**
- Prediction accuracy comparison
- TM-score plot

---

## 4.3 Selection of structures for molecular dynamics
- Explain why specific proteins were chosen
- Include:
  - gold standards
  - thermophilic candidate
  - mesophilic candidate
  - psychrophilic candidate

**Table:**
- MD simulation protein list

---

# 5. Molecular dynamics analysis of UvsX thermal adaptation

## 5.1 Simulation stability

Analyse:
- RMSD
- equilibration behaviour
- overall stability

**Figures:**
- RMSD plots

---

## 5.2 Residue flexibility

Analyse:
- RMSF profiles
- flexible regions
- differences between thermal groups

**Figures:**
- RMSF plots

---

## 5.3 Global structural properties

Analyse:

- Radius of gyration
- Compactness
- Intramolecular distances

**Figures:**
- Rg plots
- distance plots

---

## 5.4 Comparison of thermal groups

Compare:

- psychrophilic UvsX
- mesophilic UvsX
- thermophilic UvsX

Look for:
- stability differences
- flexibility differences
- structural adaptation patterns

---

# 6. Integrated identification of promising UvsX candidates

## 6.1 Combining sequence, thermal and structural evidence
- Identify candidates that satisfy multiple criteria:
  - high-confidence sequence identification
  - suitable thermal properties
  - good structural predictions
  - favourable MD behaviour

## 6.2 Prioritised candidates for future experimental testing
- Present final candidate shortlist
- Explain selection criteria

**Figures/Tables:**
- Candidate ranking table
- Final workflow figure