# Yong & Weiß et al., 2026

This repository contains all scripts and pipelines for the analyses presented in:
---
Zijuan Yong<sup>+</sup>, Josefine Friederike Weiß<sup>+</sup>, Kathleen R. Stoof-Leichsenring, Sisi Liu,Andrei Andreev, Ulrike Herzschuh<sup>*</sup> (2026): 
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

---

## Workflow Overview

1. Prepare merged sedaDNA data with taxonomic assignments and read counts per species using ['HOLI'](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025-HOLI)
2. Use [`Data_Mining/Habitat/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Data_Mining/Habitat) to generate updated species lists with habitat info.
3. Run the [`pre_genC_pipeline` Part 1](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/blob/main/genC_pipeline/pre_genC_pipeline/pre_genC_pipeline_Part1.R) to classify species based on isolation source tags and ratios and [`pre_genC_pipeline` Part 2](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/blob/main/genC_pipeline/pre_genC_pipeline/pre_genC_pipeline_Part2.R) to generate DNA weight proportions.
4. Use [`Data_Mining/DNA_C_value/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Data_Mining/DNA_C_value) to generate lists of DNA weight per cell per organism.
5. Run the [`genC_pipeline`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/genC_pipeline) to estimate group-specific DNA-based biomass and carbon.

---

## Notes

- All paths need to be adjusted to your local system.
- The project uses R and Python. See individual subfolder README files for further details.

## Author of Repository: Josefine Friederike Weiß & Zijuan Yong
## Contact: Josefine-Friederike.Weiss@awi.de
