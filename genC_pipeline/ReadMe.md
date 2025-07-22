README: genC Pipeline for Biomass Calculation
1. Overview
This R script, the genC Pipeline, is designed for the automated calculation of biomass estimates from DNA weight data across multiple sediment cores. The script is built to provide a repeatable and efficient workflow by encapsulating the processing logic into a central function and applying it iteratively to a defined list of cores.

The process includes reading raw data, proportionally distributing unclassified DNA portions, calculating 64 different biomass combinations per organism group, and finally, aggregating these results to generate statistical summaries and wide-format tables for each individual core.

2. Requirements
Required R Libraries
Ensure that the following R packages are installed before running the script:

install.packages("dplyr")
install.packages("readr")
install.packages("tidyr")
install.packages("purrr")
install.packages("tidyverse")

The script automatically loads these libraries upon execution.

Input Files
Two main data files are required to run the script:
DNA Weight Data (Processed_DNA_habitat_percentage_aqu_ter_fer_unclassfied_new.csv):
This CSV file should contain the percentage DNA portions for various taxonomic and functional groups.
It must include a Number column to identify individual cores and an Age column.
The column names for the DNA portions must end with _percentage_Start_wt.
Cell Values (values_per_cell.csv):
This CSV file contains the conversion factors (DNA_value_1 to DNA_value_8 and Biomass_value_1 to Biomass_value_8).
It requires an Organism column that maps to the organism groups defined in the script's biomass_config (e.g., "Bacteria", "Fungi", "tree").
3. Configuration and Execution
Before running the script, you need to customize three sections:

Step 1: Adjust File Paths
Update the paths to your input files and the desired output folder:

# Path to the main DNA file
DNA_weight_original <- read.csv("/PATH/TO/Processed_DNA_habitat_percentage_aqu_ter_fer_unclassfied_new.csv")

# Path to the values-per-cell file
values_per_cell <- read_delim("/PATH/TO/values_per_cell.csv", ...)

# Base path for the output results
output_base_path <- "/PATH/TO/YOUR/RESULTS_FOLDER"

Step 2: Define Cores for Processing
Customize the cores_to_process list to specify which cores should be processed. Each entry in the list requires the number (corresponding to the Number column in your CSV) and a name for naming the output folders and files.

cores_to_process <- list(
  list(number = 1, name = "Ilirney"),
  list(number = 2, name = "Salmon"),
  list(number = 3, name = "Lele"),
  list(number = 4, name = "Lama"),
  list(number = 5, name = "Btoko"),
  list(number = 6, name = "ulu")
  # Add more cores here if needed
)

Step 3: Run the Script
After the paths and the core list are configured, you can run the entire R script. It will automatically iterate through each defined core and save the results. The progress will be displayed in the console.

4. How the Script Works
Data Preprocessing: First, the two input files are loaded. Then, the unclassified DNA portions (_unclassfied_) are proportionally distributed among the main groups (aquatic, terrestrial, viruses) to update the starting weights (_Start_wt).
Main Function (process_core_biomass): This function performs the entire analysis for a single core:
It filters the data for the specified core_number.
It loops through an internal configuration list (biomass_config) that defines each organism group to be calculated.
For each group, 64 biomass combinations are calculated based on the DNA_value_ and Biomass_value_ columns.
The results are aggregated to calculate statistical metrics (mean, median, standard deviation, etc.) per age level.
Two result files are saved for the core in a dedicated subfolder.
Main Loop: An lapply loop at the end of the script calls the process_core_biomass function for each core defined in cores_to_process.
5. Output Files
For each processed core, two CSV files are created in a subfolder (named after the core_name):
[core_name]_model_biomass_sd.csv:
Contains the statistical summaries (mean, median, standard deviation, etc.) for each organism group and age level.
This file does not contain data from the "unclassfied" groups.
Format: Long format, with a Type column to distinguish between organism groups.
[core_name]_model_biomass_g_unclassfied.csv:
Contains the values for q1_biomass, mean_biomass, and median_biomass in a wide table format.
Each organism group (including the "unclassfied" groups) becomes its own set of columns.
This file is ideal for analyses where the biomass of different groups is to be compared directly.
