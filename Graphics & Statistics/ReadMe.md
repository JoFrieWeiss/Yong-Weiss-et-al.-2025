# BR_Biomass Pipeline
# Projected DNA-based Organic Carbon Burial Rate Estimation from Lake Sediments
# 
# This R pipeline estimates projected DNA-based organic carbon burial rate (BR_Biomass)
# from lake sediment cores using standardized environmental predictors and a linear mixed-effects model.
# 
# Important
# Before running this pipeline, ensure the input file 
# 'Calculate_all_information_result_unclassfied.csv' is correctly prepared and accessible.
#
# This file must contain:
# - BR_biomass (raw DNA-based burial estimate)
# - Environmental predictors: Sediment.Rate, PC1, PC2, Pann, TJul
# - Metadata: Lake, Age, Number
#
# Overview
# The pipeline includes the following main steps:
#
# 1. Load Required Libraries
#    Uses lmerTest and dplyr for data processing and mixed-effects modeling.
#
# 2. Data Input
#    - Reads main sedaDNA-environment file.
#
# 3. Preprocessing
#    - Standardizes environmental variables (z-scores) within each lake.
#    - Log10-transforms BR_biomass to improve normality.
#    - Merges standardized predictors with transformed response variable.
#
# 4. Modeling
#    - Fits a linear mixed-effects model using:
#      BR_Biomass ~ scaled_Age + scaled_Sediment.Rate + scaled_PC1 +
#                   scaled_PC2 + scaled_Pann + scaled_TJul + (1 | Lake)
#    - Lake is treated as a random intercept.
#
# 5. Output
#    - Outputs model summary using summary(fit)
#    - (Plotting and figure generation is excluded in this version)
#
# Required Input File
# - Calculate_all_information_result_unclassfied.csv
#
# Example: Running the Model
# Make sure to adjust the file path to match your environment before running the code.
#
# Authors: Zijuan Yong, Josefine Friederike Weiß
# Date: 2025
