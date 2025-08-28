################################################################################
################################################################################
############################# genC pipeline ####################################
################################################################################
# genC_pipeline.R
#
# OC DNA projected Calculation Pipeline from DNA Weight Data
#
# This script processes DNA weight data for multiple sediment cores,
# distributes unclassified DNA proportions into classified groups,
# calculates OC DNA projected values for different organism groups based on per-cell biomass parameters,
# and saves the biomass output for each core as separate CSV files.
#
# The pipeline:
# - Loads DNA weight and per-cell biomass data,
# - Normalizes and distributes unclassified DNA weights proportionally,
# - Calculates OC DNA projected values per group using predefined organism-specific parameters,
# - Iterates through all specified sediment cores,
# - Outputs processed OC DNA projected data for further analyses.
#
# Author: Zijuan Yong, Josefine Friederike Weiß
# Date: 2025-07-24


# Required Libraries
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
# library(tidypaleo) # Not used in the script
# library(ggplot2)   # Not used in the script
library(tidyverse)

# --- Define Directories ---
# 1. Path to the processed data (the output from your previous script)
processed_data_dir <- "/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/results_from_script/"

# 2. Path to the original input data
input_data_dir <- "/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/"


# --- Load Files ---

# Load the main DNA weight data from the processed data directory
tryCatch({
  DNA_weight_original  <- read.csv(paste0(processed_data_dir, "Processed_DNA_habitat_percentage_aqu_ter_fer_unclassfied_new.csv"))
  message("Main DNA file successfully loaded.")
}, error = function(e) {
  stop("ERROR: Main DNA file could not be found. Please check the path in 'processed_data_dir'.")
})


# Load the values per cell from the input directory
# NOTE: Assuming "values_per_cell.csv" is in the main directory
tryCatch({
  values_per_cell <- read_csv(paste0(input_data_dir, "values_per_cell.csv"))
  message("'values_per_cell.csv' successfully loaded.")
}, error = function(e) {
  stop("ERROR: 'values_per_cell.csv' could not be found. Please check the path in 'input_data_dir'.")
})

# Calculate the starting DNA weight and distribute the unclassified portions.
# This step is performed only once for the entire dataset.
DNA_weight <- DNA_weight_original %>%
  rowwise() %>%
  mutate(
    # Calculate sums for aquatic, terrestrial, and unclassified groups
    B_aquatic = sum(c_across(c(
      Aquatic_algae_percentage_Start_wt, Bacteria_aquatic_percentage_Start_wt,
      Viridiplantae_aquatic_percentage_Start_wt, Archaea_aquatic_percentage_Start_wt,
      Fungi_aquatic_percentage_Start_wt, Metazoa_aquatic_percentage_Start_wt
    )), na.rm = TRUE),
    
    B_terrestrial = sum(c_across(c(
      Bacteria_terrestrial_percentage_Start_wt, Viridiplantae_woody_percentage_Start_wt,
      Viridiplantae_non_woody_percentage_Start_wt, Archaea_terrestrial_percentage_Start_wt,
      Fungi_terrestrial_percentage_Start_wt, Metazoa_terrestrial_percentage_Start_wt
    )), na.rm = TRUE),
    
    B_viruses = Viruses_percentage_Start_wt,
    B_total = B_aquatic + B_terrestrial + B_viruses,
    
    B_unclassified = sum(c_across(c(
      Aquatic_algae_unclassfied_percentage_Start_wt, Bacteria_unclassfied_percentage_Start_wt,
      Viridiplantae_unclassfied_percentage_Start_wt, Archaea_unclassfied_percentage_Start_wt,
      Fungi_unclassfied_percentage_Start_wt, Metazoa_unclassfied_percentage_Start_wt,
      Viruses_unclassfied_percentage_Start_wt
    )), na.rm = TRUE),
    
    # Calculate proportions of the main groups
    P_aquatic = if_else(B_total > 0, B_aquatic / B_total, 0),
    P_terrestrial = if_else(B_total > 0, B_terrestrial / B_total, 0),
    P_viruses = if_else(B_total > 0, B_viruses / B_total, 0),
    
    # Distribute unclassified portions proportionally to the main groups
    U_aquatic = B_unclassified * P_aquatic,
    U_terrestrial = B_unclassified * P_terrestrial,
    U_viruses = B_unclassified * P_viruses,
    
    # Add the distributed portions to the respective subgroups
    # Aquatic groups
    Aquatic_algae_percentage_Start_wt = Aquatic_algae_percentage_Start_wt + if_else(B_aquatic > 0, U_aquatic * Aquatic_algae_percentage_Start_wt / B_aquatic, 0),
    Bacteria_aquatic_percentage_Start_wt = Bacteria_aquatic_percentage_Start_wt + if_else(B_aquatic > 0, U_aquatic * Bacteria_aquatic_percentage_Start_wt / B_aquatic, 0),
    Viridiplantae_aquatic_percentage_Start_wt = Viridiplantae_aquatic_percentage_Start_wt + if_else(B_aquatic > 0, U_aquatic * Viridiplantae_aquatic_percentage_Start_wt / B_aquatic, 0),
    Archaea_aquatic_percentage_Start_wt = Archaea_aquatic_percentage_Start_wt + if_else(B_aquatic > 0, U_aquatic * Archaea_aquatic_percentage_Start_wt / B_aquatic, 0),
    Fungi_aquatic_percentage_Start_wt = Fungi_aquatic_percentage_Start_wt + if_else(B_aquatic > 0, U_aquatic * Fungi_aquatic_percentage_Start_wt / B_aquatic, 0),
    Metazoa_aquatic_percentage_Start_wt = Metazoa_aquatic_percentage_Start_wt + if_else(B_aquatic > 0, U_aquatic * Metazoa_aquatic_percentage_Start_wt / B_aquatic, 0),
    
    # Terrestrial groups
    Bacteria_terrestrial_percentage_Start_wt = Bacteria_terrestrial_percentage_Start_wt + if_else(B_terrestrial > 0, U_terrestrial * Bacteria_terrestrial_percentage_Start_wt / B_terrestrial, 0),
    Viridiplantae_woody_percentage_Start_wt = Viridiplantae_woody_percentage_Start_wt + if_else(B_terrestrial > 0, U_terrestrial * Viridiplantae_woody_percentage_Start_wt / B_terrestrial, 0),
    Viridiplantae_non_woody_percentage_Start_wt = Viridiplantae_non_woody_percentage_Start_wt + if_else(B_terrestrial > 0, U_terrestrial * Viridiplantae_non_woody_percentage_Start_wt / B_terrestrial, 0),
    Archaea_terrestrial_percentage_Start_wt = Archaea_terrestrial_percentage_Start_wt + if_else(B_terrestrial > 0, U_terrestrial * Archaea_terrestrial_percentage_Start_wt / B_terrestrial, 0),
    Fungi_terrestrial_percentage_Start_wt = Fungi_terrestrial_percentage_Start_wt + if_else(B_terrestrial > 0, U_terrestrial * Fungi_terrestrial_percentage_Start_wt / B_terrestrial, 0),
    Metazoa_terrestrial_percentage_Start_wt = Metazoa_terrestrial_percentage_Start_wt + if_else(B_terrestrial > 0, U_terrestrial * Metazoa_terrestrial_percentage_Start_wt / B_terrestrial, 0),
    
    # Viruses
    Viruses_percentage_Start_wt = Viruses_percentage_Start_wt + U_viruses
  ) %>%
  ungroup()


# --- 2. Definition of the Processing Function ---

#' Processes the biomass data for a single core.
#'
#' @param core_number The 'Number' of the core to be processed.
#' @param core_name The name of the core (e.g., "Ilirney"), used for the filenames.
#' @param dna_data The complete, prepared DNA dataset.
#' @param cell_data The dataset with the values per cell.
#' @param base_output_path The base path where the resulting CSV files should be saved.
process_core_biomass <- function(core_number, core_name, dna_data, cell_data, base_output_path) {
  
  cat(paste("--- Starting processing for core:", core_name, "(Number:", core_number, ") ---\n"))
  
  # a) Filter and prepare data for the specific core
  sediment_data <- dna_data %>%
    filter(Number == core_number) %>%
    select(Age, contains("Start_wt")) %>%
    rename_with(~gsub("percentage_", "", .), .cols = everything()) %>%
    rename_with(~gsub("Start_wt", "DNA_weight", .), .cols = everything())
  
  if (nrow(sediment_data) == 0) {
    cat(paste("Warning: No data found for core", core_name, "found. Skipping.\n"))
    return(NULL)
  }
  
  age_data <- sediment_data$Age
  
  # b) Configuration for the biomass calculation
  # Defines which columns should be processed with which organisms and prefixes.
  biomass_config <- list(
    list(dna_col = "Aquatic_algae_DNA_weight", organism = "freshwater_phytoplankton", prefix = "freshwater_phytoplankton_Biomass", type = "Aquatic_algae"),
    list(dna_col = "Aquatic_algae_unclassfied_DNA_weight", organism = "freshwater_phytoplankton", prefix = "freshwater_phytoplankton_unclassfied_Biomass", type = "Aquatic_algae_unclassfied"),
    list(dna_col = "Bacteria_aquatic_DNA_weight", organism = "Bacteria", prefix = "Bacteria_aquatic_Biomass", type = "Bacteria_aquatic"),
    list(dna_col = "Bacteria_terrestrial_DNA_weight", organism = "Bacteria", prefix = "Bacteria_terrestrial_Biomass", type = "Bacteria_terrestrial"),
    list(dna_col = "Bacteria_unclassfied_DNA_weight", organism = "Bacteria", prefix = "Bacteria_unclassfied_Biomass", type = "Bacteria_unclassfied"),
    list(dna_col = "Viridiplantae_woody_DNA_weight", organism = "tree", prefix = "woody_Biomass", type = "Viridiplantae_woody"),
    list(dna_col = "Viridiplantae_non_woody_DNA_weight", organism = "nowood", prefix = "non_woody_Biomass", type = "Viridiplantae_non_woody"),
    list(dna_col = "Viridiplantae_aquatic_DNA_weight", organism = "Embryophyta", prefix = "Viridiplantae_aquatic_Biomass", type = "Viridiplantae_aquatic"),
    list(dna_col = "Viridiplantae_unclassfied_DNA_weight", organism = "Plant", prefix = "Viridiplantae_unclassfied_Biomass", type = "Viridiplantae_unclassfied"),
    list(dna_col = "Archaea_aquatic_DNA_weight", organism = "Archaea", prefix = "Archaea_aquatic_Biomass", type = "Archaea_aquatic"),
    list(dna_col = "Archaea_terrestrial_DNA_weight", organism = "Archaea", prefix = "Archaea_terrestrial_Biomass", type = "Archaea_terrestrial"),
    list(dna_col = "Archaea_unclassfied_DNA_weight", organism = "Archaea", prefix = "Archaea_unclassfied_Biomass", type = "Archaea_unclassfied"),
    list(dna_col = "Fungi_aquatic_DNA_weight", organism = "Fungi", prefix = "Fungi_aquatic_Biomass", type = "Fungi_aquatic"),
    list(dna_col = "Fungi_terrestrial_DNA_weight", organism = "Fungi", prefix = "Fungi_terrestrial_Biomass", type = "Fungi_terrestrial"),
    list(dna_col = "Fungi_unclassfied_DNA_weight", organism = "Fungi", prefix = "Fungi_unclassfied_Biomass", type = "Fungi_unclassfied"),
    list(dna_col = "Metazoa_aquatic_DNA_weight", organism = "Metazoa", prefix = "Metazoa_aquatic_Biomass", type = "Metazoa_aquatic"),
    list(dna_col = "Metazoa_terrestrial_DNA_weight", organism = "Metazoa", prefix = "Metazoa_terrestrial_Biomass", type = "Metazoa_terrestrial"),
    list(dna_col = "Metazoa_unclassfied_DNA_weight", organism = "Metazoa", prefix = "Metazoa_unclassfied_Biomass", type = "Metazoa_unclassfied"),
    list(dna_col = "Viruses_DNA_weight", organism = "Viruses", prefix = "Viruses_Biomass", type = "Viruses"),
    list(dna_col = "Viruses_unclassfied_DNA_weight", organism = "Viruses", prefix = "Viruses_unclassfied_Biomass", type = "Viruses_unclassfied")
  )
  
  # c) Loop through the configuration to calculate biomass results
  all_biomass_stats <- purrr::map_df(biomass_config, function(config) {
    
    # Check if the required DNA column exists
    if (!config$dna_col %in% names(sediment_data)) {
      cat(paste("Warning: Column", config$dna_col, "not found in the data for core", core_name, ". Skipping.\n"))
      return(NULL)
    }
    
    # Prepare data for the current organism
    dna_subset <- sediment_data %>% select(Age, !!sym(config$dna_col))
    
    cell_subset <- cell_data %>%
      filter(Organism == config$organism) %>%
      crossing(Age = age_data)
    
    merged_data <- merge(cell_subset, dna_subset, by = "Age") %>%
      mutate(!!sym(config$dna_col) := as.numeric(!!sym(config$dna_col)))
    
    # Calculate all 64 biomass combinations
    biomass_combinations <- expand.grid(dna_val = 1:8, biomass_val = 1:8)
    
    biomass_results <- merged_data
    
    for (i in 1:nrow(biomass_combinations)) {
      dna_idx <- biomass_combinations$dna_val[i]
      biomass_idx <- biomass_combinations$biomass_val[i]
      
      new_col_name <- paste0(config$prefix, "_", dna_idx, "_", biomass_idx)
      
      dna_col_name <- paste0("DNA_value_", dna_idx)
      biomass_col_name <- paste0("Biomass_value_", biomass_idx)
      
      biomass_results <- biomass_results %>%
        mutate(!!new_col_name := (!!sym(config$dna_col) / !!sym(dna_col_name)) * !!sym(biomass_col_name))
    }
    
    # d) Convert results to long format and calculate statistics
    biomass_long <- biomass_results %>%
      pivot_longer(
        cols = starts_with(config$prefix),
        names_to = "Biomass_type",
        values_to = "Biomass_value"
      )
    
    stats <- biomass_long %>%
      group_by(Age) %>%
      summarise(
        mean_biomass = mean(Biomass_value, na.rm = TRUE),
        median_biomass = median(Biomass_value, na.rm = TRUE),
        sd_biomass = sd(Biomass_value, na.rm = TRUE),
        n = n(),
        se_biomass = sd_biomass / sqrt(n),
        iqr_biomass = IQR(Biomass_value, na.rm = TRUE),
        lower_biomass = mean_biomass - (iqr_biomass / 2),
        upper_biomass = mean_biomass + (iqr_biomass / 2),
        q1_biomass = quantile(Biomass_value, 0.25, na.rm = TRUE)
      ) %>%
      mutate(Type = config$type) # Add the type for later merging
    
    return(stats)
  })
  
  # e) Save results
  # Create the specific folder for the core if it doesn't exist
  core_output_path <- file.path(base_output_path, core_name)
  if (!dir.exists(core_output_path)) {
    dir.create(core_output_path, recursive = TRUE)
  }
  
  # Save the aggregated data (without "unclassfied")
  biomass_data_filtered <- all_biomass_stats %>% filter(!grepl("unclassfied", Type))
  output_file_sd <- file.path(core_output_path, paste0(core_name, "_model_biomass_sd.csv"))
  write.csv(biomass_data_filtered, output_file_sd, row.names = FALSE)
  cat(paste("-> Results (sd) saved to:", output_file_sd, "\n"))
  
  # Save the data in wide format (with "unclassfied")
  biomass_wide <- all_biomass_stats %>%
    select(Age, Type, q1_biomass, mean_biomass, median_biomass) %>%
    arrange(Age) %>%
    pivot_wider(names_from = Type, values_from = c(q1_biomass, mean_biomass, median_biomass))
  
  output_file_wide <- file.path(core_output_path, paste0(core_name, "_model_biomass_g_unclassfied.csv"))
  write.csv(biomass_wide, output_file_wide, row.names = FALSE)
  cat(paste("-> Results (wide) saved to:", output_file_wide, "\n"))
  
  cat(paste("--- Processing for core:", core_name, "completed ---\n\n"))
  
  return(all_biomass_stats)
}


# --- 3. Main Loop for Processing All Cores ---

# Define the cores to be processed in a list
# Add or remove cores here as needed.
cores_to_process <- list(
  list(number = 1, name = "Ilirney"),
  list(number = 2, name = "Salmon"),
  list(number = 3, name = "Lele"),
  list(number = 4, name = "Lama"),
  list(number = 5, name = "Btoko"),
  list(number = 6, name = "ulu")
)

# Define the base path for the output
# NOTE: Adjust the path where the results should be saved.
output_base_path <- processed_data_dir

# Apply the function to each core in the list
# The results are stored in 'results_list' if you want to use them further in R.
results_list <- lapply(cores_to_process, function(core) {
  process_core_biomass(
    core_number = core$number,
    core_name = core$name,
    dna_data = DNA_weight,
    cell_data = values_per_cell,
    base_output_path = output_base_path
  )
})
