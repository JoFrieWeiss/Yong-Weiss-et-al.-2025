# Burial Rate and Biomass Analysis

## Overview

This repository contains the R code and processed dataset used to investigate the relationships between environmental variables, total biomass burial rate, and aquatic biomass percentage.

## Files

```text
Burial_rate_biomass_analysis/
├── README.md
├── GLMM_analysis.R
└── data/
    └── Processed_DNA_with_MC_biomass_final.csv
```

### `GLMM_analysis.R`

This script performs two generalized linear mixed-effects models (GLMMs):

1. **Total biomass burial rate**

   * Response: `BR_Total_median_biomass`
   * Gamma distribution with log link

2. **Aquatic biomass percentage**

   * Response: `Biomass_Aquatic_percentage`

Environmental predictors are standardized within Lake, and Lake is included as a random intercept in both models.

### `Processed_DNA_with_MC_biomass_final.csv`

This is the final processed dataset containing all variables required for the analyses. No additional biomass calculations are performed in `GLMM_analysis.R`.

## Required R packages

```r
install.packages(c("dplyr", "glmmTMB", "ggeffects"))
```

## Running the analysis

Place the dataset in the `data/` folder and run:

```r
source("GLMM_analysis.R")
```

This will run the complete statistical analysis.

