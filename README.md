# Taxon-resolved sources of organic carbon preserved in lake sediments using the sedaDNA GenC pipeline

This repository contains all scripts and pipelines for the analyses presented in:
---
Zijuan Yong<sup>+</sup>, Josefine Friederike Weiß<sup>+</sup>, Kathleen R. Stoof-Leichsenring, Sisi Liu, Andrei Andreev, Ulrike Herzschuh<sup>*</sup> (2026): 
Taxon-resolved sources of organic carbon preserved in lake sediments using the sedaDNA GenC pipeline


## Important prerequisites

Before using the scripts here, raw sequencing data must first be processed with the **HOLI pipeline**, which performs taxonomic assignment and data preprocessing.

You can find the HOLI pipeline and related documentation here:

- **HOLI pipeline repository:** [JoFrieWeiss/Yong-Weiss-et-al.-2025-HOLI](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025-HOLI)

Make sure to complete this step to generate the necessary input files used in the current analysis pipelines.

## Repository Structure

- [`Data_Mining/Habitat/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Data_Mining/Habitat): Scripts for scraping species-level habitat metadata (e.g., from BacDive and Microbeatlas).
- [`genC_pipeline/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline): Pipeline to estimate group-specific DNA-based biomass and carbon from sedaDNA.
- [`Pre genC pipeline/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline/Read%20length): Calculating the average read length for later processing in the genC pipeline.
- [`Post genC pipeline statistics`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Graphics%20&%20Statistics)
- [`Post genC pipeline figures`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Figures)

---

## Workflow Overview

1. Prepare merged sedaDNA data with taxonomic assignments and read counts per species using ['HOLI'](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025-HOLI)
2. Use [`Data_Mining/Habitat/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Data_Mining/Habitat) to generate updated species lists with habitat info.
3. Run the [`pre_genC_pipeline/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline/Read%20length) to extract the read length.
4. Use [`Data_Mining/DNA_C_value/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Data_Mining/DNA_C_value) to generate lists of DNA weight per cell per organism.
5. Run the [`genC_pipeline/`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/MICROSPyDER_genC_pipeline/genC%20pipeline) to estimate group-specific DNA-based biomass and carbon.
6. For further statistics please run [`Statistics`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Graphics%20&%20Statistics).
7. To reconstruct all figures go to [`Figures`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/tree/main/Figures).

---

## Notes

- All paths need to be adjusted to your local system.
- The project uses R and Python. See individual subfolder README files for further details.

## Author of Repository: Josefine Friederike Weiß & Zijuan Yong
## Contact: Josefine-Friederike.Weiss@awi.de
