# Analysis of Lake Sediment Core Data

The repository includes analyses of biomass burial rate, aquatic biomass percentage, and plant family composition using generalized linear mixed-effects models (GLMMs) and principal component analysis (PCA).

## Repository structure

```text
Yong-Weiss-et-al.-2025/
│
├── README.md
│
└── Graphics & Statistics/
    ├── GLMM_analysis.R
    ├── PCA_analysis.R
    ├── Processed_DNA_with_MC_biomass_final.csv
    └── Comibne_plant_percentage.csv
```

## Analysis scripts

### GLMM_analysis.R

This script performs generalized linear mixed-effects model analyses of:

1. Total biomass burial rate
2. Aquatic biomass percentage

Continuous predictor variables are standardized within lake, and **Lake** is included as a random intercept.

The analysis uses the following R packages:

* `dplyr`
* `glmmTMB`
* `ggeffects`

### PCA_analysis.R

This script performs principal component analysis of plant family composition.

The analysis includes:

* Hellinger transformation of plant family composition data
* Within-lake normalization
* Principal component analysis using `vegan::rda()`
* Extraction of PCA site scores
* Calculation of the variance explained by the principal components

The analysis uses the following R packages:

* `dplyr`
* `vegan`

## Input datasets

### Processed_DNA_with_MC_biomass_final.csv

This dataset contains the processed biomass and environmental variables used for the GLMM analyses, including biomass burial rate, aquatic biomass percentage, environmental predictors, and sequencing-related variables.

### Comibne_plant_percentage.csv

This dataset contains plant family composition data used for the PCA analysis.

## Reproducibility

The R scripts are designed to reproduce the statistical analyses reported in the study using the corresponding input datasets provided in this repository.

The scripts should be run from the repository root directory. Input data are located in the `Graphics & Statistics/` directory.

## R packages

The required packages can be installed using:

```r
install.packages(c(
  "dplyr",
  "glmmTMB",
  "ggeffects",
  "vegan"
))
```

## Data and code availability

The datasets and R scripts required to reproduce the statistical analyses are publicly available in this repository.

The repository contains the analysis code together with the corresponding input datasets used for the GLMM and PCA analyses.


For any questions regarding the code or the analysis, please contact Josefine Friederike Weiß at Josefine-Friederike.Weiss@awi.de
