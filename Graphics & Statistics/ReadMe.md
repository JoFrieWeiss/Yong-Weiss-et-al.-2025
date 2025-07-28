# Analysis of Lake Sediment Core Data

## Description

This project analyzes lake sediment core data to explore the relationships between organic carbon (OC) burial, biomass composition, and environmental factors over time. The primary R script (`analysis_script.R`) processes the data, performs statistical analyses (PCA, RDA, LMM), and generates a series of figures to visualize the findings.

The key analyses include:
-   Comparing Total Organic Carbon (TOC) with DNA-projected OC.
-   Visualizing the burial rates of total, aquatic, and terrestrial biomass.
-   Tracking the changing composition of taxonomic groups (e.g., bacteria, algae, plants) over the last 60,000 years.
-   Modeling the environmental drivers of OC burial rates.

---

## Prerequisites

To run this analysis, you will need:
-   **R** (version 4.0 or later recommended)
-   **RStudio** (recommended for easier project management)
-   The following R packages. You can install them by running this command in your R console:
    ```R
    install.packages(c("tidyverse", "reshape2", "cowplot", "patchwork", "scales", "vegan", "lmerTest", "effects", "mgcv", "devtools", "RColorBrewer", "ggrepel", "ggtext", "png", "grid"))
    
    # The 'corit' package needs to be installed from GitHub
    devtools::install_github("EarthSystemDiagnostics/corit")
    ```

---

## Data Files

The analysis requires the following CSV files, which should be placed in a `data` subfolder within your project directory:

1.  `data/Calculate_all_information_result_unclassfied.csv`: The main dataset containing all measurements for sediment samples, including age, TOC, DNA results, and biomass calculations.
2.  `data/Viridiplantae_PCA_combin_all.csv`: A specialized dataset containing detailed abundance data for plant families (Viridiplantae), used for the PCA and specific family time-series plots.

---

## How to Run the Script

1.  **Clone or download the repository.**
2.  **Organize your files.** Create two subdirectories in the main project folder:
    -   `data/`: Place the required `.csv` files here.
    -   `Figure/`: This folder will be created automatically by the script to store the output plots.
3.  **Update file paths.** Open the `analysis_script.R` file. At the top of the script (Section 2), modify the `data_path` and `figure_path` variables to match your local file structure. **Using relative paths is highly recommended.**
    ```R
    # Example of updated relative paths
    data_path <- "data/Calculate_all_information_result_unclassfied.csv"
    figure_path <- "Figure/"
    viridiplantae_path <- "data/Viridiplantae_PCA_combin_all.csv"
    ```
4.  **Execute the script.** Run the entire `analysis_script.R` in R or RStudio. The script will process the data and save all figures as both `.png` and `.pdf` files in the `Figure` directory.

---

## Script Structure

The R script is organized into the following sections for clarity and maintainability:

-   **Section 1: Load Libraries:** Loads all necessary R packages.
-   **Section 2: Configuration and Data Loading:** Sets file paths and loads the primary datasets. It also performs initial data cleaning and defines helper functions used throughout the script.
-   **Figure 2-6 & Supplementary Figures S3-S11:** Each figure has its own dedicated section. Within each section, the code performs the specific data manipulation and generates the plot. This modular structure makes it easy to find, understand, and re-run the code for any specific figure.

---

## Output

The script will generate a series of high-resolution plots saved in the `Figure/` directory. These include:

-   **Maps and Bar Plots:** Comparing burial rates across different lake cores.
-   **Time-Series Plots:** Showing trends in TOC, DNA, biomass, and environmental variables over time.
-   **Compositional Plots:** Visualizing the relative abundance of different taxonomic groups.
-   **Ordination Plots (PCA/RDA):** Exploring the relationships between species composition and environmental drivers.
-   **Model Diagnostic and Effect Plots (LMM/GLM):** Showing the statistical results of the linear mixed-effects models.

---

## Contact

For any questions regarding the code or the analysis, please contact [Your Name] at [your.email@example.com].
