# Analysis of Lake Sediment Core Data


## Overview

This repository contains the R scripts used for the statistical analyses in this study.

The analyses include:

* Generalized linear mixed-effects models (GLMMs) for biomass burial rate and aquatic biomass percentage
* Principal component analysis (PCA) of plant family composition
* Extraction of PCA scores for downstream analyses

## Repository structure

```text
Burial_rate_biomass_analysis/
│
├── README.md
│
├── scripts/
   ├── GLMM_analysis.R
   └── PCA_analysis.R

```

## Analyses

### GLMM analysis

`GLMM_analysis.R` contains two generalized linear mixed-effects models:

* **Total biomass burial rate**

  * Response: `BR_Total_median_biomass`
  * Gamma distribution with log link

* **Aquatic biomass percentage**

  * Response: `Biomass_Aquatic_percentage`

Environmental variables are standardized within Lake, and Lake is included as a random intercept.

### PCA analysis

`PCA_analysis.R` performs PCA on plant family composition data.

The analysis includes:

1. Hellinger transformation
2. Standardization within Lake
3. PCA using `vegan::rda()`
4. Extraction of PCA site scores
5. Calculation of explained variance

## Data availability

The datasets used in these analyses are not included in this repository.

The scripts require the corresponding input datasets to run the analyses. Data availability is subject to the data-sharing conditions of the study.

## R packages

The analyses require the following R packages:

```r
install.packages(c(
  "dplyr",
  "glmmTMB",
  "ggeffects",
  "vegan"
))
```

## Running the analyses

The scripts are located in the `scripts/` directory.

Run the scripts from the repository root:

```r
source("scripts/GLMM_analysis.R")
source("scripts/PCA_analysis.R")
```

The scripts reproduce the statistical analyses used in this study.


## Contact

For any questions regarding the code or the analysis, please contact Josefine Friederike Weiß at Josefine-Friederike.Weiss@awi.de
