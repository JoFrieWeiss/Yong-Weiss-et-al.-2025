# Figures 2–7

R scripts used to generate Figures 2–7 for the analysis of DNA-projected biomass and organic carbon dynamics.

## Figure overview

| Figure       | Description                                              | Script       |
| ------------ | -------------------------------------------------------- | ------------ |
| **Figure 2** | OC DNA-projected burial rate                             | `Figure_2.R` |
| **Figure 3** | DNA, biomass, TOC, burial rate, temperature, and PC1     | `Figure_3.R` |
| **Figure 4** | TOC vs. DNA-projected biomass                            | `Figure_4.R` |
| **Figure 5** | Aquatic vs. terrestrial biomass and Botryococcus         | `Figure_5.R` |
| **Figure 6** | DNA weight and DNA-projected biomass composition         | `Figure_6.R` |
| **Figure 7** | OC DNA-projected biomass percentage over the last 30 kyr | `Figure_7.R` |

## Data

All figures use the same processed dataset:

```text
data/Processed_DNA_with_MC_biomass_final.csv
data/Pollendata/Algae pecentage csv.csv
```

## Requirements

The scripts were written in R and use packages including:

```text
dplyr
tidyr
ggplot2
scales
patchwork
cowplot
zoo
stringr
ggthemes
readr
```

## Running the scripts

Run the scripts from the project root directory:

```r
source("Figure_2.R")
source("Figure_3.R")
source("Figure_4.R")
source("Figure_5.R")
source("Figure_6.R")
source("Figure_7.R")
```

Output figures are saved in:

```text
figures
```
## Notes
The datasets and R scripts required to reproduce Figures 2–7 are publicly available in this repository.
For questions regarding the code or analyses, please contact Zijuan Yong at Zijuanyong@163.com.
