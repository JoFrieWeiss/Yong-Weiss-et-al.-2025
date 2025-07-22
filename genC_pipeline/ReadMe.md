# genC Pipeline

Biomass Estimation from Sedimentary Ancient DNA (sedaDNA)  
This R pipeline estimates group-specific biomass (e.g. aquatic, terrestrial, viruses) from sedimentary ancient DNA (sedaDNA) using DNA weight proportions and conversion values per cell.

## Overview

The pipeline includes the following main steps:

1. **Load Required Libraries**  
   Uses `tidyverse`, `readr`, `dplyr`, and `purrr` for data processing.

2. **Data Input**  
   - Main sedaDNA file: contains relative DNA weights for different organism groups across sediment cores.
   - Cell-specific conversion table: average DNA weight and biomass per cell for various organism types.

3. **Preprocessing**  
   - Calculates the total DNA per group (aquatic, terrestrial, viruses).
   - Distributes unclassified DNA proportionally to the classified groups.
   - Updates DNA weight per subgroup accordingly.

4. **Per-Core Biomass Calculation**  
   For each sediment core:
   - Filters the core-specific DNA data.
   - Matches it with corresponding cell values by organism.
   - Applies all 64 possible DNA-to-biomass conversion factor combinations.
   - Computes summary statistics (mean, median, standard deviation, interquartile range, etc.) by depth/age.

5. **Output**  
   - Saves two output CSVs per core:
     - A detailed file with biomass statistics per group.
     - A wide-format summary for downstream analysis or plotting.

## Usage

Make sure to adapt all file paths in the script to your local or server environment.

### Required Input Files
- `Processed_DNA_habitat_percentage_aqu_ter_fer_unclassfied_new.csv`  
  Relative DNA weight per taxonomic group per sample.

- `values_per_cell.csv`  
  Contains average DNA weight and biomass per cell for various organism types. Make sure decimal separators are correct (e.g. `,` if European format).

### Example: Run for All Cores

You can define a list of sediment cores like this:

```r
cores_to_process <- list(
  list(number = 1, name = "Ilirney"),
  list(number = 2, name = "Salmon"),
  list(number = 3, name = "Lele"),
  list(number = 4, name = "Lama"),
  list(number = 5, name = "Btoko"),
  list(number = 6, name = "ulu")
)
