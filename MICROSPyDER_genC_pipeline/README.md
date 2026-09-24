# MICROSPyDER & genC Pipeline
**Microbial Community Source Tracking & Stochastic OC Reconstruction**

This repository contains the integrated pipeline for ecological habitat classification (**MICROSPyDER**) and stochastic OC reconstruction (**genC**). The model is specifically designed to handle environmental shotgun metagenomics data, accounting for laboratory biases and biological variability.

---

## 1. MICROSPyDER: Ecological Source Tracking
The MICROSPyDER module classifies the ecological origin of microbial taxa by calculating the **Habitat Affinity Score (HAS)**. This score contrasts the frequency of a taxon $i$ in terrestrial vs. aquatic environments based on global MicrobeAtlas data.

### HAS Formula
The HAS score ($T_i$) for a specific taxon $i$ is defined as:

$$T_i = \frac{\sum_{j \in \mathcal{H}_{terr}} V_{i,j} - \sum_{k \in \mathcal{H}_{aqua}} V_{i,k}}{\sum_{j \in \mathcal{H}_{terr}} V_{i,j} + \sum_{k \in \mathcal{H}_{aqua}} V_{i,k}}$$

**Key Parameters:**
* $V_{i,j}$: Frequency of taxon $i$ in terrestrial habitat $j$.
* $V_{i,k}$: Frequency of taxon $i$ in aquatic habitat $k$.
* ℋ<sub>terr</sub> and ℋ<sub>aqua</sub>: Mathematical sets representing the collection of terrestrial and aquatic habitat descriptors.
        
**Key Parameters:**
* $T_i \in [-1, 1]$
* **Classification:** Specialists are identified at $|T_i| \geq 0.8$, while generalists cluster around $0.0$.

---

## 2. genC: OC Reconstruction
The genC pipeline reconstructs the physical total OC $\hat{B}_i$ using a stochastic **Monte Carlo Integration** ($N=10,000$). It transforms sequencing read counts into organic carbon mass per gram of sediment.

### OC Formula
The reconstructed OC for taxon $i$ is calculated through the global equation:

$$\hat{B}_{i} = \mathbb{M}_{s=1 \dots N} \left[ \frac{ \left( \text{MW} \cdot [DNA]_{p} \cdot \bar{L} \cdot 10^{-3} \right) \cdot \text{CDF} }{ \text{CCF} \cdot (M_{raw} \cdot (1 - \theta)) } \cdot R_i \cdot \frac{ C_{i,s} \cdot S_i \cdot E }{ D_{i,s} } \right]$$

**Variables and Units:**
* $\hat{B}_{i}$: Median reconstructed OC ($g_{C} \cdot g_{sed}^{-1}$)
* $\bar{L}$: Arithmetic mean read length of the sample (bp)
* $[DNA]_{p}$: DNA concentration in the sequencing pool ($ng \cdot \mu l^{-1}$)
* $R_i$: Relative read proportion (%)
* $C_{i,s}$: Simulated carbon mass per cell via Beta-Priors ($g_{C} \cdot cell^{-1}$)
* $D_{i,s}$: Simulated DNA mass per cell via Beta-Priors ($g_{DNA} \cdot cell^{-1}$)
* $E$: Extraction efficiency factor (fixed at 3.5)
* $S_i$: Structural correction factor (e.g., 3.0 for Bacteria, 15.0 for Woody Plants)
* $MW$: Molecular weight of nucleotides
* $CDF / CCF$: Lab dilution and concentration factors
* $M_{raw} \cdot (1 - \theta)$: Dry sediment mass calculation

---

## 3. Repository Structure
* [`1. MICROSPyDER using HAS`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline/1_v2_TAS_GenC_pipeline.R): HAS score calculation.
* [`2. Calculating DNA weight`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline/2_v2_2_DNAweight_GenC.R): Merging sequencing metadata with taxonomic data to calculate DNA weight per taxonomic group.
* [`3. Calculating OC`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline/3_v2_3_GenC.R): Stochastic MC simulation and linear TOC-comparison plots.
* [`Data used`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/Data): Contains reference files (e.g., `values_per_cell_correct.csv`).

## 4. Usage
1. Adjust the `base_dir` in all R-scripts to match your local path.
2. Run the scripts sequentially (Part 1 $\rightarrow$ Part 2 $\rightarrow$ Part 3).
3. Ensure R packages `dplyr`, `ggplot2`, `patchwork`, and `sensitivity` are installed.

## 5. Sensitivity Analysis & Robustness
To evaluate the reliability of the reconstructed OC $\hat{B}_i$, the pipeline includes a global sensitivity analysis (Monte Carlo based). This identifies which input parameters (e.g., DNA concentration, C-value, or extraction efficiency) contribute most to the overall uncertainty.

### Mathematical Approach: OAT & Sobol Indices
We use a combination of **One-at-a-Time (OAT)** variation and **Variance-based Sensitivity Analysis**. 

The total variance of the output $V(Y)$ is decomposed into the contributions of individual input factors:

$$V(Y) = \sum_{i} V_i + \sum_{i < j} V_{ij} + \dots + V_{12\dots k}$$

* **First-order Sensitivity Index ($S_i$):** Measures the direct contribution of input $X_i$ to the output variance:
    $$S_i = \frac{V(E[Y | X_i])}{V(Y)}$$
* **Total-effect Index ($S_{Ti}$):** Accounts for the main effect of $X_i$ plus all its higher-order interactions with other variables.

---
**Authors:** Josefine Friederike Weiß (2026)  
**Project:** Nature Communications (Preparation)
