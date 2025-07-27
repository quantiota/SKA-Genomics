# SKA Genome: Information Geometry of the Human Genome

**Exploring the hidden information architecture of the genome using Structured Knowledge Accumulation (SKA).**



## Overview

**SKA-Genome** applies the SKA learning framework to DNA sequences, generating information-theoretic “entropy maps” of the genome.  
Unlike traditional genomics tools, SKA analyzes the *information architecture* of sequences in a forward-only, model-free fashion—highlighting novel, functionally significant, or complex regions without the need for prior annotation or motif databases.



## Why SKA?

- **Model-free:** No need for predefined motifs, gene models, or references.
- **Adaptive:** SKA “learns” directly from the data stream—mapping entropy and knowledge along the sequence.
- **Universal:** Works for DNA, RNA, protein, epigenomic marks, or even time-series derived from single-cell or tumor data.
- **Unsupervised discovery:** Flags informational anomalies, boundaries, and regime changes invisible to variant-centric or alignment-based methods.


## Getting Started

1. Download a FASTA file (e.g., human chromosome)
2. Run the SKA script:
   ```bash
   python src/ska_genome_analysis.py --input data/Homo_sapiens.GRCh38.dna.chromosome.21.fa --output figures/chr21_entropy.png

## Quick Start (Demo)

For rapid testing or educational purposes, we recommend beginning with the classic *E. coli* K-12 reference genome:

- **Small size** (4.6 Mb) for fast prototyping
- **Universal compatibility** with all genomics tools
- **Extensively annotated** for comparison and validation

Download:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
```




## Example Applications

- **Genome annotation:** Detect genes, regulatory elements, repetitive regions, and structural boundaries.
- **Cancer genomics:** Compare normal and tumor tissue; track clonal evolution or chromosomal instability.
- **Evolutionary genomics:** Identify conserved vs. rapidly evolving regions; discover species-specific “information signatures.”
- **Personalized medicine:** Explore patient-specific entropy landscapes; highlight novel or rare genomic events.
- **Single-cell & multi-omics:** Map heterogeneity and transitions in complex cell populations.



## Folder Structure

- `data/` — Genome or sequence files (FASTA/CSV)
- `src/` — SKA algorithm and utilities
- `figures/` — Output plots and example results
- `README.md` — Project documentation and instructions



## Scientific Impact

SKA-Genome enables:

- Unsupervised mapping of the genome’s information landscape at any scale (nucleotide, gene, chromosome, whole genome)

- Discovery of boundaries and functional elements even in poorly annotated or novel genomes

- Comparative information profiling across samples, tissues, disease states, or species

- Insights into cancer, evolution, and genome organization beyond classical variant-centric methods

## Foundation Theory

The SKA-Genome project is based on the theoretical framework and mathematical foundations of Structured Knowledge Accumulation (SKA) developed by QUANTIOTA.

For in-depth theory, proofs, and reference papers, please see our open-access arXiv repository:

[QUANTIOTA / Arxiv — Foundation Theory](https://github.com/quantiota/Arxiv)

This resource covers the core principles that power SKA’s universal, model-free learning across genomics, physics, finance, and more.



## Citation

If you use this work, please cite:  
> **Bouarfa Mahi, “SKA-Genome: Information Geometry of the Human Genome” (2025), GitHub.**






