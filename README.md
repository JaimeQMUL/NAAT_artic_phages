# UvsX Ortholog Identification - Cold Adaptation Research

## 🎯 Project Overview

This repository contains the complete bioinformatics pipeline for identifying **cold-adapted UvsX orthologs** across various microbial species. The project focuses on:

- Curating comprehensive UvsX and RecA protein databases
- Developing Hidden Markov Models (HMMs) for domain detection
- Analyzing temperature profiles for cold-adaptation signatures
- Constructing phylogenetic trees with domain annotations
- Generating comprehensive visualizations for publication

## 📊 Latest Achievements

**Current Status:** ✅ **All major analysis stages complete**

**Most Recent Progress (Commit: `d927287`):**
- ✅ Successfully identified temperature profiles of **novel UvsX homologs**
- ✅ Completed **T4 UvsX analysis**
- ✅ Sandbox analysis identified additional homologs
- ✅ Novel homolog discovery with temperature data

## 🔬 Research Pipeline

### 1. Data Collection & Curation ✅

**Completed Tasks:**
- Pulled all UvsX-like and RecA-like protein orthologs from NCBI and UniProt
- Created comprehensive curated database
- Extracted sequences using Pfam domain annotations
- Removed redundant sequences
- Scripted recreation of lab-engineered mutants (KAN-5, KAN-7, KAN-9)
- All mutants from literature recreated programmatically

**Data Sources:**
- NCBI Protein Database
- UniProt Knowledgebase
- Pfam Domain Database
- Literature-curated mutant sequences

### 2. Sequence Processing ✅

**Tools:**
- `cleaning_data.py` - Data preprocessing and quality control
- `filtering_data.py` - Sequence filtering and redundancy removal
- Modular data wrangling utilities
- Directory structure organization

**Pipeline:**
```
Raw Sequences → Quality Filtering → Redundancy Removal → Curated Database
```

### 3. HMM Analysis ✅

**Implementation:**
- Created HMM models for UvsX and RecA domain detection
- Integrated HMMER3 for sensitive domain searches
- Generated optimized HMM profiles
- Temperature profile optimization using `seq2topt`
- Novel homolog identification with HMM validation

**Key Features:**
- High-sensitivity domain detection
- Temperature-specific profile optimization
- Cold-adaptation signature identification
- Novel sequence discovery

### 4. Multiple Sequence Alignment ✅

**Tools:**
- `alignment.py` - MSA generation
- `phylotree.py` - Phylogenetic tree construction
- Tree visualization and domain mapping

**Outputs:**
- Aligned sequence matrices
- Phylogenetic trees with domain annotations
- Tree-based evolutionary analysis

### 5. Statistical Analysis ✅

**Methods:**
- Shannon entropy calculations
- Family/genus/species distribution analysis
- Sequence length distributions
- Statistical validation of cold-adaptation signatures

**Metrics:**
- Entropy-based adaptation scoring
- Statistical significance testing
- Distribution analysis across taxonomic ranks

### 6. Complete Visualizations ✅

**Generated Figures:**
- Family counts (`family_counts.png`)
- Genus counts (`genus_counts.png`)
- Species counts (`species_counts.png`)
- Ranks summary (`ranks_summary.png`)
- Sequence length distributions (`sequence_lengths.png`)
- Tree visualizations with domain mapping (`mapped_only_blue_tree.png`)
- Temperature profile plots
- Complete figure suite in `uvsx/plots/Jaime_ESPrint_Protein_Visualisation/`

## 🔧 Technical Architecture

### Project Structure

```
Dissertation/
├── src/                    # Core bioinformatics tools
│   ├── biotools/          # FASTA, MSA utilities
│   │   └── fasta_tools.py
│   │   └── msa.py
│   ├── extracting_references/  # Database creation from UniProt
│   │   └── creating_curated_database.py
│   ├── helpers/           # Data wrangling, file creation
│   │   ├── data_wrangling.py
│   │   └── directory_creation.py
│   ├── pfam/              # Pfam domain tools
│   │   └── pfam.py
│   └── visualisations/    # Plot generation tools
│       ├── fasta_summary_plots.py
│       ├── phylogenetic_trees.py
│       └── shannon_entropy.py
├── seq2topt/              # Temperature optimization models
│   ├── __init__.py
│   ├── functions.py       # Core analytical functions
│   ├── model.py           # ML models for cold-adaptation
│   └── seq2topt.py        # Main temperature analysis script
├── uvsx/                  # UvsX-specific analysis
│   ├── data/              # Raw and processed sequence data
│   └── plots/             # Generated figures and visualizations
│       ├── Jaime_ESPrint_Protein_Visualisation/
│       │   └── Figures/
│       ├── family_counts.png
│       ├── genus_counts.png
│       ├── species_counts.png
│       ├── ranks_summary.png
│       └── mapped_only_blue_tree.png
├── notebooks/             # Jupyter analysis notebooks
│   ├── data_wrangling.ipynb
│   ├── domains.ipynb
│   ├── hmm_practice.ipynb
│   ├── temperatures.ipynb
│   └── notebook.ipynb
├── hmm/                   # HMM models
│   └── *.hmm
├── docs/                  # Documentation
│   └── documentation/
│       ├── cleaning_database.md
│       ├── creating_curated_database.md
│       ├── domains.md
│       └── hmms.md
├── alignment.py          # MSA generation
├── cleaning_data.py      # Data preprocessing
├── filtering_data.py     # Quality filtering
├── phylotree.py          # Phylogenetic tree building
├── uvsx_database.py      # Database queries
├── summary_visualisations.py
└── environment.yml       # Conda environment
```

### Dependencies

```yaml
# environment.yml
name: dissertation-env
channels:
  - conda-forge
dependencies:
  - python=3.10
  - biopython
  - hhmmer
  - numpy
  - pandas
  - matplotlib
  - seaborn
  - scikit-learn
```

## 🧬 Key Research Outcomes

### Completed Milestones

✅ **Database Creation:** Complete curated UvsX ortholog database  
✅ **HMM Models:** Validated models for cold-adaptation detection  
✅ **Phylogenetics:** Multiple trees with domain annotations  
✅ **Temperature Analysis:** Cold-adaptation profile calculations  
✅ **Novel Homologs:** Identified new UvsX sequences with temperature data  
✅ **Visualizations:** Complete figure suite for publication  
✅ **Reproducible Pipeline:** Modular, documented, and testable code  

### Research Findings

- Successfully identified temperature profiles of novel UvsX homologs
- Completed comprehensive T4 UvsX analysis
- Discovered additional homologs through sandboxed testing
- Established cold-adaptation signatures across multiple species
- Validated HMM models for temperature-specific domain detection

## 🚀 How to Use

### Setup

```bash
# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate dissertation-env
```

### Running Analysis

```bash
# Clean and process data
python cleaning_data.py

# Filter sequences
python filtering_data.py

# Run HMM analysis
python seq2topt/seq2topt.py

# Generate visualizations
python src/visualisations/phylogenetic_trees.py
```

### Generating Plots

```bash
# Run all visualizations
python summary_visualisations.py
```

## 📊 Project Status

| Component | Status | Notes |
|-----------|--------|-------|
| Data Collection | ✅ Complete | Curated database ready |
| HMM Development | ✅ Complete | Models validated |
| Phylogenetic Analysis | ✅ Complete | Trees generated |
| Temperature Analysis | ✅ Complete | Cold-adaptation signatures identified |
| Visualizations | ✅ Complete | All figures generated |
| Novel Homolog Discovery | ✅ Complete | Novel sequences identified |
| Documentation | ✅ Complete | Methodology documented |

## 📝 Next Steps

- Analyze results for cold-adaptation signatures
- Refine analysis pipeline based on findings
- Prepare dissertation chapters
- Publication preparation
- Further species expansion (if needed)

## 📚 References

- HMMER3 documentation
- Pfam domain database
- UniProt Knowledgebase
- NCBI Protein Database
- Literature on UvsX cold-adaptation

## 🔗 Contact

For questions or collaboration: [Your Contact Information]

## 📄 License

[Add appropriate license information]

---

*This project is part of a Master's dissertation on identifying cold-adapted UvsX orthologs across microbial species.*
