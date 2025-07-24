# pre_genC_pipeline.R
#
# Pre-Processing Pipeline for DNA and Habitat Data Integration and Calculation of Derived Metrics
#
# This script loads complex DNA information and habitat percentage data,
# merges them, converts relevant columns to numeric, calculates DNA- and sediment-related derived metrics,
# and creates new columns representing weighted DNA contributions by habitat percentages.
#
# Author: Josefine Friederike Weiß, Zijuan Yong
# Date: 2025-07-24

library(readr)
library(dplyr)

# --- Functions --------------------------------------------------------

#' Load CSV with minimal name repair and optional suppression of column type messages
#'
#' @param path Character, file path to CSV file
#' @param suppress_col_types Logical, whether to suppress printing column type info (default FALSE)
#' @return Data frame read from CSV
load_csv <- function(path, suppress_col_types = FALSE) {
  read_csv(path, name_repair = "minimal", show_col_types = !suppress_col_types)
}

#' Convert specified columns of a data frame to numeric safely
#'
#' @param df Data frame
#' @param cols Character vector of column names to convert
#' @return Data frame with specified columns converted to numeric (if they exist)
convert_to_numeric <- function(df, cols) {
  existing_cols <- intersect(cols, names(df))
  df[existing_cols] <- lapply(df[existing_cols], function(x) as.numeric(as.character(x)))
  return(df)
}

#' Compute derived DNA and sediment metrics based on input data frame columns
#'
#' @param df Data frame containing required columns
#' @return Data frame with new columns appended
compute_derived_metrics <- function(df) {
  df %>%
    mutate(
      Total_concentration = `Total concentration` * 1000, # Convert to consistent unit
      Final_DNA_in_Pool = (`Molecular weight nucleotide` * `Dilution concentration before sequencing 2` * `Average read length of sample`) * 0.001,
      Cumulated_Dilution_Factor = Total_concentration / `Dilution concentration before sequencing 2`,
      Cumulated_Concentration_Factor = `Library concentration` / `Concentration before Genejet`,
      Starting_DNA_ng = (Final_DNA_in_Pool * Cumulated_Dilution_Factor) / Cumulated_Concentration_Factor,
      Starting_DNA_g = Starting_DNA_ng * 1e-9,
      Dry_weight = `DNA` * (1 - `Water content` * 0.01),
      Starting_DNA_wt = Starting_DNA_g / Dry_weight,
      ARBULK = `Sediment Rate` * `Dry Bulk Density`,
      BRTOC = ARBULK * TOC * 0.01
    )
}

#' Create weighted DNA columns by multiplying habitat percentage columns with Starting_DNA_wt
#'
#' @param df Data frame
#' @param pct_cols Character vector of habitat percentage columns
#' @return Data frame with new columns appended: "{original_column}_Start_wt"
create_start_wt_columns <- function(df, pct_cols) {
  for (col in pct_cols) {
    new_col <- paste0(col, "_Start_wt")
    if (col %in% names(df)) {
      df[[new_col]] <- df[[col]] * df$Starting_DNA_wt
    } else {
      warning(paste("Column missing:", col))
      df[[new_col]] <- NA_real_
    }
  }
  new_cols <- paste0(pct_cols, "_Start_wt")
  df <- convert_to_numeric(df, new_cols)
  return(df)
}

# --- Main Execution ---------------------------------------------------

# File paths - adjust these paths before running
dna_info_path <- "/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/Complex DNA information.csv"
habitat_data_path <- "/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/combine_habitat_percentage_aqu_ter_fre_unclassfied_new.csv"
output_path <- "/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/results_from_script/Processed_DNA_habitat_percentage_aqu_ter_fer_unclassfied_new_fine.csv"

# Load datasets
dna_info <- load_csv(dna_info_path)
habitat_data <- load_csv(habitat_data_path, suppress_col_types = TRUE)

# Merge on "Number" and "Age", keep only one "Core" column
combined_data <- merge(dna_info, habitat_data, by = c("Number", "Age")) %>%
  select(-Core.y) %>%
  rename(Core = Core.x)

# Define columns expected to hold percentages and numerics
percentage_columns <- c(
  "Bacteria_aquatic_percentage", "Bacteria_terrestrial_percentage", "Bacteria_unclassfied_percentage",
  "Archaea_aquatic_percentage", "Archaea_terrestrial_percentage", "Archaea_unclassfied_percentage",
  "Fungi_aquatic_percentage", "Fungi_terrestrial_percentage", "Fungi_unclassfied_percentage",
  "Viridiplantae_aquatic_percentage", "Viridiplantae_woody_percentage", "Viridiplantae_non_woody_percentage",
  "Viridiplantae_unclassfied_percentage",
  "Metazoa_aquatic_percentage", "Metazoa_terrestrial_percentage", "Metazoa_unclassfied_percentage",
  "Viruses_percentage", "Viruses_unclassfied_percentage",
  "Aquatic_algae_percentage", "Aquatic_algae_unclassfied_percentage"
)

numeric_columns <- c(
  "Total concentration", "Molecular weight nucleotide", "Dilution concentration before sequencing 2",
  "Average read length of sample", "Library concentration", "Concentration before Genejet",
  "DNA", "Water content", "Sediment Rate", "Dry Bulk Density", "TOC"
)

# Convert relevant columns to numeric to avoid issues in calculations
combined_data <- convert_to_numeric(combined_data, numeric_columns)
combined_data <- convert_to_numeric(combined_data, percentage_columns)

# Compute DNA and sediment derived metrics
combined_data <- compute_derived_metrics(combined_data)

# Create weighted DNA columns based on habitat percentages
combined_data <- create_start_wt_columns(combined_data, percentage_columns)

# Save processed results
write.csv(combined_data, output_path, row.names = FALSE)

message("Processing completed successfully. Output saved to: ", output_path)
