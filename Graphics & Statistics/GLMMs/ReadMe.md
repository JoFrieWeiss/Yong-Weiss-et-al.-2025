# Projected DNA-based Organic Carbon Burial Rate Estimation from Lake Sediments

This R pipeline estimates the projected DNA-based organic carbon burial rate (`BR_Biomass`) from lake sediment cores using standardized environmental predictors and a linear mixed-effects model.

---

## Prerequisites

Before running this pipeline, ensure that the input file [`Calculate_all_information_result_unclassfied.csv`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/blob/main/Graphics%20%26%20Statistics/Calculate_all_information_result_unclassfied.csv)  is correctly prepared and accessible.

This file must contain the following columns:
- **BR_biomass**: raw DNA-based burial rate estimate  
- **Environmental predictors**: `Sediment.Rate`, `PC1`, `PC2`, `Pann`, `TJul`  
- **Metadata**: `Lake`, `Age`, `Number`  

---

## Pipeline Overview

The pipeline consists of the following main steps:

### 1. Load Required Libraries
Loads the `lmerTest` and `dplyr` packages for statistical modeling and data processing.

### 2. Data Input
Reads the input file containing sedaDNA and environmental data.

### 3. Preprocessing
- Standardizes (z-scores) the environmental predictors within each lake.
- Log10-transforms `BR_biomass` to improve normality.
- Merges standardized predictors with the transformed response variable.

### 4. Modeling
Fits a linear mixed-effects model with the following formula:

```R
BR_Biomass ~ scaled_Age + scaled_Sediment.Rate + scaled_PC1 +
             scaled_PC2 + scaled_Pann + scaled_TJul + (1 | Lake)

