# Analysis of Lake Sediment Core Data

## Description

This project analyzes lake sediment core data to explore the relationships between organic carbon (OC) burial, biomass composition, and environmental factors over time. The project contains two main R scripts:

1.  **[`Plots.R`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/blob/main/Graphics%20%26%20Statistics/Plots.R)**: This script focuses on data processing, exploration, and the generation of all primary and supplementary figures for publication.
2.  **[`PCA_RDA.R`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/blob/main/Graphics%20%26%20Statistics/PCA_RDA.R)**: This script is dedicated to performing focused multivariate statistical analyses, including Redundancy Analysis (RDA) and Principal Component Analysis (PCA).

---

## Prerequisites

To run this analysis, you will need:
-   **R** (version 4.0 or later recommended)
-   **RStudio** (recommended for easier project management)
-   The following R packages. You can install all of them by running this command in your R console:
    ```R
    # This command installs all packages needed for both scripts
    install.packages(c(
        "tidyverse", "reshape2", "cowplot", "patchwork", "scales", 
        "vegan", "lmerTest", "effects", "mgcv", "devtools", 
        "RColorBrewer", "ggrepel", "ggtext", "png", "grid"
    ))
    
    # The 'corit' package needs to be installed from GitHub
    devtools::install_github("EarthSystemDiagnostics/corit")
    ```

---

## Data Files

The analysis requires the following CSV files, which should be placed in a `data` subfolder within your project directory:

1. [`Calculate_all_information_result_unclassfied.csv`](https://github.com/JoFrieWeiss/Yong-Weiss-et-al.-2025/blob/main/Graphics%20%26%20Statistics/Calculate_all_information_result_unclassfied.csv) : The main dataset containing all measurements for sediment samples, including age, TOC, DNA results, and biomass calculations. Used by both scripts.
2.  [`Viridiplantae_PCA_combin_all.csv`]: A specialized dataset containing detailed abundance data for plant families (Viridiplantae), used for the PCA in both scripts.

---

## Project Scripts and Usage

This project is organized around two primary R scripts.

### 1. Main Visualization Script (`analysis_script.R`)

This is the primary script for creating the visual outputs of the project.

-   **Purpose**: To process the raw data and generate all figures (e.g., Fig 2-6, Supplementary Figs S3-S11).
-   **Usage**:
    1.  Ensure all required data files are in the `data` folder.
    2.  Update the file paths at the top of the script.
    3.  Run the script. It will generate and save all plots to the `Figure` folder.
-   **Output**: A series of high-resolution plots saved as both `.png` and `.pdf` files.

### 2. Statistical Analysis Script (`multivariate_analysis.R`)

This script is dedicated to the core statistical modeling and does not produce saved plots.

-   **Purpose**: To perform focused multivariate statistical analyses on the dataset. The script is structured to run RDA and PCA, and is set up with packages for GLMM (Generalized Linear Mixed-Effects Models).
-   **Key Methods**:
    -   **Redundancy Analysis (RDA)**: Performed on the full taxonomic dataset. The species data (biomass percentages) is **fourth-root transformed** to stabilize variance before analysis.
    -   **Principal Component Analysis (PCA)**: Performed on the Viridiplantae (plant) dataset. The species data undergoes a **Hellinger transformation** followed by a **group-wise normalization** for each lake.
-   **Usage**: Run the script to perform the calculations. The results, such as model summaries, will be printed directly to the R console.
-   **Output**: Console output containing summaries of the RDA and PCA models (e.g., explained variance, component scores).

---

## Output Summary

-   The **`analysis_script.R`** generates a comprehensive set of visual outputs (plots) in the `Figure/` directory.
-   The **`multivariate_analysis.R`** provides statistical model results directly in the R console for interpretation and further analysis.

---

## Contact

For any questions regarding the code or the analysis, please contact Josefine Friederike Weiß at Josefine-Friederike.Weiss@awi.de
