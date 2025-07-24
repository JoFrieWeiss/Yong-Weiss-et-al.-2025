# Yong & Weiß et al., 2025

This repository contains all scripts and pipelines for the analyses presented in:
---
Zijuan Yong<sup>+</sup>, Josefine Friederike Weiß<sup>+</sup>, Kathleen R. Stoof-Leichsenring, Sisi Liu, Ulrike Herzschuh<sup>*</sup> (2025): 
Quantitative contribution of major organism groups to lake organic carbon burial since the last glacial inferred from sedimentary ancient DNA


## Important prerequisites

Before using the scripts here, raw sequencing data must first be processed with the **HOLI pipeline**, which performs taxonomic assignment and data preprocessing.

You can find the HOLI pipeline and related documentation here:

- **HOLI pipeline repository:** [JoFrieWeiss/Yong-Weiss-et-al.-2025-HOLI](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025-HOLI)

Make sure to complete this step to generate the necessary input files used in the current analysis pipelines.

## Repository Structure

- [`Data_Mining/Habitat/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Data_Mining/Habitat): Scripts for scraping species-level habitat metadata (e.g., from BacDive).
- [`genC_pipeline/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/genC_pipeline): Pipeline to estimate group-specific DNA-based biomass and carbon from sedaDNA.
  - [`pre_genC_pipeline/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/genC_pipeline/pre_genC_pipeline): Preprocessing to calculate DNA weights and prepare habitat proportions.
    - `Habitat_Classification/`: Scripts to classify taxa into habitat types using isolation metadata and ratios.
    - `Taxonomic_Assignment/`: Scripts to process sedaDNA taxonomic data (e.g., from HOLI).

---

## Workflow Overview

1. Prepare merged sedaDNA data with taxonomic assignments and read counts per species.
2. Use `Data_Mining/Habitat/BacDive.py` to generate updated species lists with habitat info.
3. Use `Habitat_Classification/` to classify species based on isolation source tags and ratios.
4. Run the `pre_genC_pipeline` to generate DNA weight proportions.
5. Run the `genC_pipeline` to estimate group-specific DNA-based biomass and carbon.

---

## Notes

- All paths need to be adjusted to your local system.
- The project uses R and Python. See individual subfolder README files for further details.
