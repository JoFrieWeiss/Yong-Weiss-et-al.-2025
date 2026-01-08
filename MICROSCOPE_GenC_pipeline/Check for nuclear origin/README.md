# Bioinformatic Pipeline: Taxon-Specific Extraction & Re-mapping

This pipeline is designed for high-precision taxonomic verification. It consists of two interconnected SLURM scripts: **Step 1** filters broad metaDMG results into specific taxonomic subsets, and **Step 2** performs targeted re-alignment of those specific reads against chosen reference genomes.

[Image of bioinformatic workflow showing taxonomic read filtering from metaDMG results followed by Bowtie2 mapping]

---

## 🛠 SETUP & PREREQUISITES

### 1. Cluster Modules
The scripts require the following tools to be available in your environment:
- `bowtie2/2.5.1`
- `samtools/1.20`

### 2. Metadata File (`sample2org2ref.tsv`)
Step 2 requires a tab-separated metadata file located in your working directory. This file tells the script which reference genome to use for each organism found in your samples:

```text
sample	organism	ref_genome	ref_size
Sample_A	Salicaceae	GCF_00123	1500
Sample_A	Dryas	GCF_00456	2000
