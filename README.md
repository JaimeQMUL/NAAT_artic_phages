# UvsX Ortholog Identification and Cold-Adaptation Research

## Project Overview

This repository contains the bioinformatics and computational analysis pipeline developed for a Master's dissertation investigating whether **cold-adapted UvsX orthologs** can be identified across microbial species.

The project follows four main stages:


1. Database Curation
        ↓
2. Hidden Markov Model Analysis
        ↓
3. In-Silico Thermal Screening
        ↓
4. Molecular Dynamics


The overall aim is to identify UvsX homologs, evaluate their predicted temperature adaptation, and investigate whether selected candidates show structural and dynamic properties consistent with cold adaptation.

---

## Research Pipeline

### 1. Database Curation

`1-DatabaseCuration/`

Collection and processing of UvsX-like sequences from NCBI, UniProt and InterPro. Sequences are cleaned, filtered and organised into curated databases, with experimentally characterised proteins maintained as gold standards.

### 2. Hidden Markov Models

`2-HiddenMarkovModels/`

Development and application of HMM profiles to identify and validate UvsX/RecA-like proteins and relevant domains. HMM benchmarking and related investigations are also contained in `Analysis/`.

### 3. In-Silico Screening

`3-InSilicoScreening/`

Prediction and scoring of candidate protein temperature profiles to identify potential cold-adapted UvsX orthologs. Top candidates are selected for downstream structural analysis.

### 4. Molecular Dynamics

`4-MolecularDynamics/`

Structural and molecular-dynamics analysis of selected candidates and gold-standard proteins using predicted and crystal structures and GROMACS. Analyses include measures such as RMSD, RMSF and radius of gyration.

---

## Repository Structure


dissertation/
│
├── 1-DatabaseCuration/       # Sequence collection and database curation
├── 2-HiddenMarkovModels/     # HMM development and analysis
├── 3-InSilicoScreening/      # Thermal prediction and candidate screening
├── 4-MolecularDynamics/      # Molecular dynamics and structural analysis
│
├── Analysis/                 # Supporting, benchmarking and exploratory analyses
├── docs/                     # Documentation
├── src/                      # Shared source code
├── writeup/                  # Dissertation writing
│
├── environment.yml           # Conda environment
└── README.md


---

## Data Flow


NCBI / UniProt / InterPro
          ↓
   Database Curation
          ↓
      HMM Analysis
          ↓
  Thermal Screening
          ↓
 Candidate Selection
          ↓
 Molecular Dynamics
          ↓
Structural / Dynamic Analysis


---

## Current Status

The main computational pipeline has been established, including:

- Curated UvsX sequence databases
- HMM-based homolog identification
- Temperature-based candidate screening
- Candidate ranking and selection
- Protein structure prediction
- Molecular dynamics workflows
- Structural and evolutionary analyses
- Supporting HMM and structure-prediction benchmarking

The remaining work focuses primarily on interpretation of the results, refinement of analyses, and dissertation/publication preparation.

---

## Installation

Create the Conda environment using:

```bash
conda env create -f environment.yml
conda activate cold
```

Individual pipeline stages can be run from their respective numbered directories. See the scripts and documentation within each directory for stage-specific instructions.

---

## Research Question

> **Can naturally occurring UvsX orthologs with signatures of cold adaptation be identified computationally and characterised through sequence, structural and molecular-dynamics analyses?**

This project aims to connect **sequence conservation, predicted temperature adaptation, evolutionary relationships, and structural dynamics** to identify promising cold-adapted UvsX candidates.