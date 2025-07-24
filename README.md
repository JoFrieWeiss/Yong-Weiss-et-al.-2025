# Yong-Weiss-et-al.-2025

This repository contains all scripts and pipelines for the analyses presented in:

# ** Zijuan Yong1,4+, Josefine Friederike Weiß1+, Kathleen R. Stoof-Leichsenring1, Sisi Liu1, Ulrike Herzschuh123* (2025): 
# Quantitative contribution of major organism groups to lake organic carbon burial since the last glacial inferred from sedimentary ancient DNA**

---

## Repository Structure

- `Data_Mining/Habitat/`: Scripts for scraping species-level habitat metadata (e.g., from BacDive).
- `genC_pipeline/`: Pipeline to estimate group-specific DNA-based biomass and carbon from sedaDNA.
  - `pre_genC_pipeline/`: Preprocessing to calculate DNA weights and prepare habitat proportions.
- `Habitat_Classification/`: Scripts to classify taxa into habitat types using isolation metadata and ratios.
- `Taxonomic_Assignment/`: Scripts to process sedaDNA taxonomic data (e.g., from HOLI).
- `Metadata/`: Supplementary species metadata (e.g., BacDive URLs, isolation source info).

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
