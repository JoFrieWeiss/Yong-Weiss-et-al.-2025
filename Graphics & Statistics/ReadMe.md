Project Overview
This project analyzes the projected OC DNA-projected burial rate across multiple lake sediment cores using environmental predictors. The main objective is to quantify how sedimentological and climatic variables influence BR_Biomass while accounting for spatial variation among lakes.
Important
Before running the BR_Biomass pipeline, make sure your input data file is prepared and located correctly:
Calculate_all_information_result_unclassfied.csv

This input must include:

DNA burial estimates (BR_biomass)

Environmental variables (e.g., Sediment Rate, PC1/PC2, Pann, TJul)

Sample metadata (e.g., Lake, Age)
Overview
The pipeline includes the following key steps:

1. Load Required Libraries
Uses lmerTest, dplyr, and readr for model fitting and data processing.

2. Data Input
Reads the main input file Calculate_all_information_result_unclassfied.csv.

3. Preprocessing
Selects relevant environmental and DNA variables.

Standardizes environmental predictors (z-scores) within each lake.

Transforms BR_biomass using log10 to improve normality.

Merges datasets by sample ID and age.

4. Modeling
Fits a linear mixed-effects model:
BR_Biomass ~ scaled_Age + scaled_Sediment.Rate + scaled_PC1 + scaled_PC2 + scaled_Pann + scaled_TJul + (1 | Lake)

Lake is treated as a random effect to account for spatial structure.

5. Output
Provides a summary of model coefficients and significance (via summary()).

Designed for statistical inference and downstream plotting (not included in this script version).

Required Input Files
Calculate_all_information_result_unclassfied.csv
Contains all required sedaDNA burial and environmental variables.

Usage Notes
Make sure all file paths in the script are adapted to your system.

Install all required packages beforehand:
install.packages(c("lmerTest", "dplyr"))


install.packages(c("lmerTest", "dplyr"))
Authors:
Zijuan Yong, Josefine Friederike Weiß
Date: 2025
