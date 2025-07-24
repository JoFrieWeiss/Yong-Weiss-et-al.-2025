# Load necessary libraries
library(tidyverse)
library(data.table)
library(readxl)
library(readr)
library(tidyr)
library(stringr)

# ---- Load data ----
# Read updated bacteria/archaea species list (already scraped habitat information using the data mining script)
#'/Volumes/projects/biodiv/user/ziyong/1_ilirney Lake/species_list/updated_Archaea_Ilirney_list.csv' 
#'#'/Volumes/projects/biodiv/user/ziyong/1_ilirney Lake/species_list/updated_Bacteria_Ilirney_list.csv'
#'
updated_Bacteria_list_Ilirney <- read_csv("/Volumes/projects/biodiv/user/ziyong/1_ilirney Lake/species_list/updated_Bacteria_Ilirney_list.csv") 

# Filter merged dataset for Bacteria at species level and age > 0
# Load merged data from HOLI taxonomic assignment
# Create a temporary environment
temp_env <- new.env()

# Load the .RData file into the temporary environment
load("/Volumes/projects/biodiv/user/ziyong/A_shotgun_data/ilirney_lca.v1_MergedData.RData", envir = temp_env)

# List objects in the file
ls(temp_env)

# Assume the object inside is called "MergedData" — rename it when assigning
Ilirney_MergedData <- temp_env$MergedData

# filter for bacteria/archaea 

Bacteria_Ilirney <- Ilirney_MergedData %>%
  filter(superkingdom == "Bacteria", rank == "species", years > 0)

# ---- Calculate habitat ratios ----
updated_Bacteria_list_Ilirney <- updated_Bacteria_list_Ilirney %>%
  mutate(
    Total_sample_number = aquatic_samples + soil_samples + animal_samples + plant_samples,
    ratio_terrestrial = soil_samples / Total_sample_number,
    ratio_aquatic = aquatic_samples / Total_sample_number,
    ratio_animals = animal_samples / Total_sample_number,
    ratio_plants = plant_samples / Total_sample_number
  )

# Merge read data with habitat metadata
combined_Bacteria_Ilirney <- merge(Bacteria_Ilirney, updated_Bacteria_list_Ilirney, by = "species")

# ---- Habitat decision function ----
habitat_decision_making <- function(data) {
  data %>%
    mutate(
      habitat_1 = case_when(
        # Freshwater environments → aquatic
        str_detect(isolation_categories, regex("#(freshwater|lake|river|pond)", ignore_case = TRUE)) ~ "aquatic",
        
        # Marine → marine
        str_detect(isolation_categories, regex("#marine", ignore_case = TRUE)) ~ "marine",
        
        # Terrestrial environments → terrestrial
        str_detect(isolation_categories, regex("#(terrestrial|soil|sediment|mud|dust)", ignore_case = TRUE)) ~ "terrestrial",
        
        TRUE ~ NA_character_
      ),
      
      habitat = case_when(
        # Use already classified
        habitat_1 %in% c("aquatic", "terrestrial", "marine") ~ habitat_1,
        
        # Mixed tags involving freshwater + terrestrial → aquatic
        str_detect(isolation_categories, regex("#(freshwater|lake|river|pond)", ignore_case = TRUE)) &
          str_detect(isolation_categories, regex("#terrestrial", ignore_case = TRUE)) ~ "aquatic",
        
        # Ratio-based classification
        ratio_terrestrial > 0.2 & ratio_aquatic > 0.2 ~ "aquatic/terrestrial",
        ratio_terrestrial > 0.2 ~ "terrestrial",
        ratio_aquatic > 0.2 ~ "aquatic",
        ratio_animals > 0.1 ~ "animal_associated",
        ratio_plants > 0.1 ~ "plant_associated",
        
        TRUE ~ NA_character_
      )
    )
}

# ---- Relative abundance functions ----
process_Bacteria_data_habitat <- function(data) {
  data %>%
    group_by(habitat, years) %>%
    summarise(sumcount = sum(taxonReads), .groups = "drop") %>%
    group_by(years) %>%
    mutate(rel_abund = sumcount / sum(sumcount) * 100) %>%
    pivot_wider(names_from = habitat, values_from = rel_abund, values_fill = 0) %>%
    pivot_longer(-years, names_to = "habitat", values_to = "rel_abund")
}

process_Bacteria_data_species <- function(data) {
  data %>%
    group_by(species, years) %>%
    summarise(sumcount = sum(taxonReads), .groups = "drop") %>%
    group_by(years) %>%
    mutate(rel_abund = sumcount / sum(sumcount) * 100) %>%
    pivot_wider(names_from = species, values_from = rel_abund, values_fill = 0) %>%
    pivot_longer(-years, names_to = "species", values_to = "rel_abund")
}

# ---- Run classification and export ----
combined_Bacteria_Ilirney_habitat <- habitat_decision_making(combined_Bacteria_Ilirney)

combined_Bacteria_Ilirney_species_long <- process_Bacteria_data_species(combined_Bacteria_Ilirney_habitat)

# Merge with habitat info
small_combined_Bacteria_Ilirney_habitat <- combined_Bacteria_Ilirney_habitat %>%
  select(species, habitat, habitat_1, isolation_categories, bacdive_url)

combined_Bacteria_Ilirney_species_long_merged <- merge(
  combined_Bacteria_Ilirney_species_long,
  small_combined_Bacteria_Ilirney_habitat,
  by = "species"
) %>%
  distinct()

# Export results
write_csv(combined_Bacteria_Ilirney_species_long_merged,
          "/Volumes/projects/biodiv/user/ziyong/1_Ilirney Lake/species_list/results_from_script/combined_Bacteria_Ilirney_species_long_merged_dedupe.csv")

write_csv(combined_Bacteria_Ilirney_habitat,
          "/Volumes/projects/biodiv/user/ziyong/1_Ilirney Lake/species_list/results_from_script/combined_Bacteria_Ilirney_habitat.csv")
