# SKA-Genome: Information Geometry of the Human Genome

**Exploring the hidden information architecture of the genome using Structured Knowledge Accumulation (SKA).**



## Overview

This repository applies the SKA learning framework to genomic sequences—revealing regions of order, disorder, and evolutionary significance through a streaming, model-free information approach.

- **Input:** Raw genomic DNA sequence (FASTA or string)
- **Output:** Entropy/knowledge profiles highlighting genes, repeats, regulatory motifs, and transitions



## Why SKA?

Unlike traditional models, SKA requires no training, annotation, or prior knowledge—just the sequence itself.  
It provides:
- A forward-only, one-pass “learnability landscape” across the genome
- Adaptive, unsupervised discovery of structure and novelty
- A tool for comparing genomes, identifying functionally important regions, and mapping evolutionary change



## Getting Started

1. **Clone this repository**
2. **Place or download a FASTA file or DNA sequence in `data/`**
3. **Run `notebooks/SKA_genome_analysis.ipynb` to generate and plot entropy/knowledge profiles**



## Example Applications

- Identifying gene-rich or gene-poor regions
- Detecting repetitive elements, conserved motifs, and regulatory sequences
- Comparing entropy landscapes across species or individuals



## Folder Structure

- `data/` — Genome or sequence files (FASTA/CSV)
- `notebooks/` — Analysis notebooks
- `src/` — SKA algorithm and utilities
- `figures/` — Output plots and example results



## Citation

If you use this work, please cite:  
> **Bouarfa Mahi, “SKA-Genome: Information Geometry of the Human Genome” (2025), GitHub.**



**Contributions and datasets welcome!**


