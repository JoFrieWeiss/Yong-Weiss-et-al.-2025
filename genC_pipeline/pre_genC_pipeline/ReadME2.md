# Complete sedaDNA Habitat Analysis Pipeline

## Overview

This repository contains a single, unified R script (`master_pipeline.R`) that performs a complete habitat analysis for sedimentary ancient DNA (sedaDNA) data. The pipeline processes data from multiple sediment cores, classifies taxa from all kingdoms (Prokaryotes and Eukaryotes) into habitat categories, and calculates their relative abundances.

The entire workflow is automated within this single script, progressing from raw input files to a final, analysis-ready summary table.

## Features

* **All-in-One Execution**: Runs the entire workflow from start to finish with a single command.
* **Multi-Core Processing**: Automatically processes data for multiple sediment cores in a loop.
* **Robust Path Management**: Uses a centralized configuration section with lookup tables to handle inconsistent file names.
* **Specialized Classification**: Implements distinct, tailored logic for:
    * **Prokaryotes** (Bacteria, Archaea) using data-mined habitat information.
    * **Eukaryotes** (Fungi, Viridiplantae, Metazoa) using a curated database and a set of refinement rules.
* **Reproducible Workflow**: Creates intermediate, classified files for prokaryotes before the final aggregation step, ensuring the workflow is transparent and reproducible.
* **Single Tidy Output**: Generates one comprehensive CSV file with all results, ready for downstream analysis and plotting in tools like `ggplot2`.

---

## Required Files

Before running the script, ensure you have the following input files:

1.  **Core Data Files (`.RData`)**: One file for each sediment core, containing a `MergedData` object from the HOLI taxonomic assignment pipeline.
2.  **Prokaryote Habitat Lists (`.csv`)**: Separate CSV files for Bacteria and Archaea for each core. These files must contain pre-scraped habitat information (e.g., from BacDive), including columns like `species`, `isolation_categories`, `aquatic_samples`, `soil_samples`, etc.
3.  **Eukaryote Habitat Database (`.xlsx`)**: A single Excel file mapping Eukaryote species names to their primary habitat.
4.  **Plant Family Database (`.csv`)**: A single CSV file mapping plant families to a type (e.g., `tree`, `no wood`) to classify terrestrial plants.

---

## Configuration & Usage

Follow these steps to run the pipeline:

### 1. Configure Paths in `SECTION 1`

Open the R script. The **only section** you need to edit is `SECTION 1: CONFIGURATION`.

Update the named vectors (`core_rdata_paths`, `bacteria_habitat_paths`, `archaea_habitat_paths`) to provide the full, correct file path for each of your input files. The names you provide (e.g., `"Ilirney"`, `"Salmon"`) will be used as identifiers in the final output.

```R
# Example Configuration in SECTION 1.3

core_rdata_paths <- c(
  Ilirney = "/path/to/your/ilirney_lca.v1_MergedData.RData",
  Salmon  = "/path/to/your/salmon_lca.v2_MergedData.RData"
  # ... and so on for all cores
)

bacteria_habitat_paths <- c(
  Ilirney = "/path/to/your/updated_Bacteria_Ilirney_list.csv",
  Salmon  = "/path/to/your/updated_Bacteria_Salmon_list.csv"
  # ... and so on for all cores
)
