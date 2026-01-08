# Bioinformatic Pipeline: From metaDMG to Genomic Visualization

This pipeline provides a complete workflow for extracting taxon-specific reads, re-mapping them against high-quality reference genomes, and visualizing the distribution across Nuclear, Mitochondrial, and Plastid DNA.

---

## 📋 Table of Contents
1. [Overview](#-overview)
2. [Prerequisites](#-prerequisites)
3. [Step 1: Taxon Extraction (SLURM)](#-step-1-taxon-extraction-slurm)
4. [Step 2: Read Re-mapping (SLURM)](#-step-2-read-re-mapping-slurm)
5. [Step 3: R Visualization (Local)](#-step-3-r-visualization-local)

---

## 🔍 Overview
The pipeline solves the problem of verifying metaDMG hits by re-aligning identified reads to specific references. This helps to confirm the presence of specific taxa and see where the reads originate (e.g., organellar vs. nuclear).

---

## 🛠 Prerequisites

### Cluster (Linux/SLURM)
- **Modules:** `bowtie2/2.5.1`, `samtools/1.20`
- **Metadata:** A file named `sample2org2ref.tsv` with columns: `sample`, `organism`, `ref_genome`, `ref_size`.

### Local Machine (R)
- **Libraries:** `tidyverse`, `fs`
- **Data:** Sequence Reports (`.tsv`) downloaded from the NCBI Assembly database for each reference genome used.

---
