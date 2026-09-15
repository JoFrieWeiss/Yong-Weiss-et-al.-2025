# Analysis of Lake Sediment Core Data

## Overview

This repository contains the R scripts used for the statistical analyses presented in this study.

The analyses include:

* Generalized linear mixed-effects models (GLMMs) for biomass burial rate and aquatic biomass percentage
* Principal component analysis (PCA) of plant family composition
* Extraction of PCA scores for further analyses

## Repository structure

```text
Burial_rate_biomass_analysis/
│
├── README.md
│
└── Graphics & Statistics/
    ├── GLMM_analysis.R
    └── PCA_analysis.R
```

## Analysis scripts

### GLMM_analysis.R

This script performs two GLMM analyses:

1. Total biomass burial rate
2. Aquatic biomass percentage

Environmental variables are standardized within Lake, with Lake included as a random intercept.

### PCA_analysis.R

This script performs PCA on plant family composition data, including:

* Hellinger transformation
* Standardization within Lake
* PCA using `vegan::rda()`
* Extraction of PCA site scores
* Calculation of explained variance

## Data availability

The datasets used for these analyses are not included in this repository.

The scripts require the corresponding input datasets to run the analyses. Data availability is subject to the data-sharing conditions of the study.

## Required R packages

```r
install.packages(c(
  "dplyr",
  "glmmTMB",
  "ggeffects",
  "vegan"
))
```

## Reproducibility

The scripts in `Graphics & Statistics/` contain the statistical analysis procedures used in this study.

## Contact

For any questions regarding the code or the analysis, please contact Josefine Friederike Weiß at Josefine-Friederike.Weiss@awi.de
