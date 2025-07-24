# pre_genC_pipeline

This R pipeline automates the processing and calculation of DNA and sediment data from multiple raw input files. It merges datasets, converts relevant columns to numeric types, computes various derived metrics, and generates a fully processed output table with key DNA and habitat-related values.

---

## Overview

The pipeline performs the following key steps:

- Reads DNA concentration and sequencing metadata from CSV files
- Reads habitat composition percentages for different taxa groups
- Merges these datasets on common keys (`Number` and `Age`)
- Converts specified columns to numeric data types to ensure proper calculations
- Calculates derived variables such as starting DNA amounts, dry weight normalization, sediment bulk density, and carbon content
- Creates additional columns representing DNA amounts weighted by habitat percentage contributions
- Saves the processed and augmented dataset as a CSV file for downstream analysis

---

## Functions in the Script

- **`load_csv()`**: Loads a CSV file with minimal column name repair and optional suppression of column type messages.
- **`convert_to_numeric()`**: Safely converts selected columns to numeric, handling missing columns gracefully.
- **`compute_derived_metrics()`**: Calculates key derived DNA and sediment metrics using relevant input columns.
- **`create_start_wt_columns()`**: Generates new columns multiplying habitat percentage columns by the normalized starting DNA weight.

---

## Input Data

- **DNA information file**: Contains raw measurements such as DNA concentration, molecular weight, dilution factors, sequencing read lengths, library concentrations, water content, sedimentation rates, bulk density, and total organic carbon.
- **Habitat percentage file**: Contains percentage values for aquatic, terrestrial, and unclassified fractions across multiple taxa groups including Bacteria, Archaea, Fungi, Viridiplantae, Metazoa, Viruses, and aquatic algae.

---

## Usage Instructions

1. **Set file paths**  
   Edit the script variables `dna_info_path`, `habitat_data_path`, and `output_path` to point to your local data files and desired output location.

2. **Install dependencies**  
   Ensure you have the following R packages installed:
   ```r
   install.packages(c("readr", "dplyr"))

