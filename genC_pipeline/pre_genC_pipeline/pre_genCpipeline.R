#----------------------------------------------------------------#
# 1. SETUP: Load libraries and global data
#----------------------------------------------------------------#

# Load required libraries
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(purrr) # Used for iterating and managing lists

# --- Load global helper data only once ---
message("Loading global helper data...")
final_results_eukaryota_excel <- read_excel("/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/final_results_eukaryota_excel.xlsx") %>%
  rename(species = Species)

plant_family_types <- read_csv("/Volumes/projects/biodiv/user/ziyong/A_Eukaryota_Habitat/Plant_Family_Type.csv")

dna_information <- read_csv("/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/Complex DNA information.csv", name_repair = "minimal")

# --- Define constants for the pipeline ---
# **FIXED**: Added specific 'file_name' for each core.
CORES_CONFIG <- list(
  Ilirney = list(id = 1, path = "/Volumes/projects/biodiv/user/ziyong/1_ilirney Lake/", 
                 file_name = "ilirney_lca.v1_MergedData.RData", filter_age = 4.4),
  Salmon  = list(id = 2, path = "/Volumes/projects/biodiv/user/ziyong/2_Salmon Lake/", 
                 file_name = "salmon_lca.v2_MergedData.RData", filter_age = c(0, 3.735)), # <-- Corrected filename
  Lele    = list(id = 3, path = "/Volumes/projects/biodiv/user/ziyong/3_Lele Lake/", 
                 file_name = "lele_lca.v1_MergedData.RData", filter_age = 4.7),
  Lama    = list(id = 4, path = "/Volumes/projects/biodiv/user/ziyong/4_Lama Lake/", 
                 file_name = "Lama_lca_MergedData.RData", filter_age = NULL),
  Btoko   = list(id = 5, path = "/Volumes/projects/biodiv/user/ziyong/5_Btoko Lake/", 
                 file_name = "btoko_lca.v1_MergedData.RData", filter_age = NULL),
  Ulu     = list(id = 6, path = "/Volumes/projects/biodiv/user/ziyong/6_Ulu Lake/", 
                 file_name = "ulu_lca.v1_MergedData.RData", filter_age = 35.81)
)


# New age values for the "Salmon" core
SALMON_NEW_AGES <- c(0, 0.163, 0.602, 1.094, 2.515, 3.622, 3.735, 3.851, 6.876, 7.489, 8.48, 10.175, 9.928, 10.282, 10.985, 11.919, 12.56, 13.056, 13.847, 14.332, 15.569, 18.467, 20.656, 22.359, 23.006, 24.868, 26.866, 28.725, 29.491)

# Base path for outputs
OUTPUT_BASE_PATH <- "/Volumes/projects/biodiv/user/ziyong/A_aquatic_terrestrial/"

#----------------------------------------------------------------#
# 2. MAIN PROCESSING FUNCTION
#----------------------------------------------------------------#
process_core <- function(core_name, config) {
  
  message(paste0("\nProcessing Core: ", core_name, " ----------------\n"))
  
  # --- 2.1 Load and Prepare Core Data ---
  # **FIXED**: Now uses the specific file_name from the config list.
  file_to_load <- file.path(config$path, config$file_name)
  message(paste("Loading file:", file_to_load))
  load(file_to_load)
  
  core_data_filtered <- MergedData %>% filter(ka > 0)
  
  total_counts <- core_data_filtered %>%
    group_by(ka) %>%
    summarise(sumcount_all = sum(taxonReads, na.rm = TRUE), .groups = "drop")
  
  super_na_counts <- core_data_filtered %>%
    filter(is.na(superkingdom)) %>%
    group_by(ka) %>%
    summarise(super_na_sumcount = sum(taxonReads, na.rm = TRUE), .groups = "drop")
  
  base_counts <- total_counts %>%
    left_join(super_na_counts, by = "ka") %>%
    mutate(super_na_sumcount = ifelse(is.na(super_na_sumcount), 0, super_na_sumcount),
           sumcount_exna = sumcount_all - super_na_sumcount)
  
  # --- 2.2 Eukaryota Pre-processing ---
  eukaryota_classified <- core_data_filtered %>%
    filter(superkingdom == "Eukaryota", rank == "species") %>%
    inner_join(final_results_eukaryota_excel, by = "species") %>%
    mutate(
      new_habitat = case_when(
        phylum == "Streptophyta" & !family %in% c("Butomaceae", "Characeae", "Closteriaceae", "Desmidiaceae", "Elatinaceae", "Equisetaceae", "Hydrocharitaceae", "Juncaceae", "Nymphaeaceae", "Lentibulariaceae", "Potamogetonaceae", "Typhaceae", "Juncaginaceae") ~ "terrestrial",
        order == "Lepidoptera" ~ "terrestrial",
        family == "Potamogetonaceae" ~ "freshwater",
        species == "Allium schoenoprasum" ~ "terrestrial",
        species %in% c("Equisetum hyemale", "Luzula luzuloides", "Sparganium erectum", "Triglochin maritima") ~ "freshwater/terrestrial",
        species == "Penaeus chinensis" | order == "Octopoda" ~ "marine",
        family == "Cordylidae" ~ "terrestrial",
        genus == "Cherax" | species %in% c("Boleophthalmus pectinirostris", "Salvelinus fontinalis") ~ "freshwater",
        TRUE ~ Habitat
      )
    )
  
  # --- 2.3 Process by Taxonomic Group ---
  
  calculate_eukaryota_percentages <- function(kingdom_name) {
    df <- eukaryota_classified %>% filter(kingdom == kingdom_name)
    
    kingdom_total <- df %>% group_by(ka) %>% summarise(!!paste0(kingdom_name, "_sumcount_all") := sum(taxonReads, na.rm = TRUE))
    
    habitat_sum <- df %>%
      mutate(habitat_group = case_when(
        new_habitat %in% c("freshwater", "marine", "marine/freshwater") ~ "aquatic",
        new_habitat == "terrestrial" | new_habitat == "freshwater/terrestrial" ~ "terrestrial",
        TRUE ~ "varied"
      )) %>%
      group_by(ka, habitat_group) %>%
      summarise(sumcount = sum(taxonReads, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = habitat_group, values_from = sumcount, values_fill = 0)
    
    combined <- base_counts %>%
      left_join(kingdom_total, by = "ka") %>%
      left_join(habitat_sum, by = "ka") %>%
      mutate(across(everything(), ~replace_na(.x, 0)))
    
    if (kingdom_name == "Viridiplantae") {
      plant_type_sum <- df %>%
        filter(new_habitat == "terrestrial") %>%
        left_join(plant_family_types, by = "family") %>%
        group_by(ka, Type) %>%
        summarise(sumcount = sum(taxonReads, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Type, values_from = sumcount, values_fill = 0, names_repair = "minimal") %>%
        rename(woody = any_of("tree"), non_woody = any_of("no wood"))
      
      combined <- combined %>% 
        left_join(plant_type_sum, by = "ka") %>%
        mutate(across(c(woody, non_woody), ~replace_na(.x, 0))) %>%
        mutate(
          total_classified_reads = aquatic + woody + non_woody,
          !!paste0(kingdom_name, "_woody_percentage") := ifelse(sumcount_all > 0, woody / sumcount_all, 0),
          !!paste0(kingdom_name, "_non_woody_percentage") := ifelse(sumcount_all > 0, non_woody / sumcount_all, 0),
          !!paste0(kingdom_name, "_terrestrial_percentage") := 0
        )
    } else {
      combined <- combined %>%
        mutate(
          total_classified_reads = aquatic + terrestrial,
          !!paste0(kingdom_name, "_terrestrial_percentage") := ifelse(sumcount_all > 0, terrestrial / sumcount_all, 0),
          !!paste0(kingdom_name, "_woody_percentage") := 0,
          !!paste0(kingdom_name, "_non_woody_percentage") := 0
        )
    }
    
    final_result <- combined %>%
      mutate(
        !!paste0(kingdom_name, "_na_count") := (!!sym(paste0(kingdom_name, "_sumcount_all")) / sumcount_exna) * super_na_sumcount,
        !!paste0(kingdom_name, "_unclassified_count") := (!!sym(paste0(kingdom_name, "_sumcount_all")) - total_classified_reads) + !!sym(paste0(kingdom_name, "_na_count")),
        !!paste0(kingdom_name, "_aquatic_percentage") := ifelse(sumcount_all > 0, aquatic / sumcount_all, 0),
        !!paste0(kingdom_name, "_unclassified_percentage") := ifelse(sumcount_all > 0, (!!sym(paste0(kingdom_name, "_unclassified_count")) / sumcount_all), 0),
        !!paste0(kingdom_name, "_percentage") := ifelse(sumcount_all > 0, (!!sym(paste0(kingdom_name, "_sumcount_all")) / sumcount_all), 0)
      ) %>%
      select(ka, starts_with(kingdom_name))
    
    return(final_result)
  }
  
  calculate_prokaryote_percentages <- function(superkingdom_name) {
    habitat_file <- file.path(config$path, "species_list", paste0("combined_", superkingdom_name, "_", core_name, "_habitat.csv"))
    
    if (!file.exists(habitat_file)) {
      message(paste("Warning: Habitat file not found for", superkingdom_name, "in", core_name, "- skipping."))
      return(NULL)
    }
    
    habitat_data <- read_csv(habitat_file, show_col_types = FALSE) %>%
      filter(ka > 0) %>%
      mutate(new_habitat = coalesce(habitat_1, habitat)) %>%
      group_by(ka, tax_id) %>%
      slice_max(order_by = taxonReads, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    superkingdom_total <- core_data_filtered %>%
      filter(superkingdom == superkingdom_name) %>%
      group_by(ka) %>%
      summarise(!!paste0(superkingdom_name, "_sumcount_all") := sum(taxonReads, na.rm = TRUE))
    
    habitat_sum <- habitat_data %>%
      group_by(ka, new_habitat) %>%
      summarise(sumcount = sum(taxonReads, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = new_habitat, values_from = sumcount, values_fill = 0)
    
    base_counts %>%
      left_join(superkingdom_total, by = "ka") %>%
      left_join(habitat_sum, by = "ka") %>%
      mutate(across(everything(), ~replace_na(.x, 0))) %>%
      mutate(
        ter_aqu_sumcount = terrestrial + aquatic,
        !!paste0(superkingdom_name, "_na_count") := (!!sym(paste0(superkingdom_name, "_sumcount_all")) / sumcount_exna) * super_na_sumcount,
        !!paste0(superkingdom_name, "_unclassified_count") := (!!sym(paste0(superkingdom_name, "_sumcount_all")) - ter_aqu_sumcount) + !!sym(paste0(superkingdom_name, "_na_count")),
        !!paste0(superkingdom_name, "_aquatic_percentage") := ifelse(sumcount_all > 0, aquatic / sumcount_all, 0),
        !!paste0(superkingdom_name, "_terrestrial_percentage") := ifelse(sumcount_all > 0, terrestrial / sumcount_all, 0),
        !!paste0(superkingdom_name, "_unclassified_percentage") := ifelse(sumcount_all > 0, (!!sym(paste0(superkingdom_name, "_unclassified_count")) / sumcount_all), 0),
        !!paste0(superkingdom_name, "_percentage") := ifelse(sumcount_all > 0, !!sym(paste0(superkingdom_name, "_sumcount_all")) / sumcount_all, 0)
      ) %>%
      select(ka, starts_with(superkingdom_name))
  }
  
  viruses_results <- core_data_filtered %>%
    filter(superkingdom == "Viruses") %>%
    group_by(ka) %>%
    summarise(Viruses_sumcount_all = sum(taxonReads, na.rm = TRUE)) %>%
    right_join(base_counts, by = "ka") %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    mutate(
      Viruses_na_count = (Viruses_sumcount_all / sumcount_exna) * super_na_sumcount,
      Viruses_percentage = Viruses_sumcount_all / sumcount_all,
      Viruses_unclassified_percentage = Viruses_na_count / sumcount_all
    ) %>%
    select(ka, starts_with("Viruses"))
  
  aquatic_algae_results <- core_data_filtered %>%
    filter(superkingdom == "Eukaryota", is.na(kingdom)) %>%
    group_by(ka) %>%
    summarise(Aquatic_algae_sumcount_all = sum(taxonReads, na.rm = TRUE)) %>%
    right_join(base_counts, by = "ka") %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    mutate(
      Aquatic_algae_na_count = (Aquatic_algae_sumcount_all / sumcount_exna) * super_na_sumcount,
      Aquatic_algae_percentage = Aquatic_algae_sumcount_all / sumcount_all,
      Aquatic_algae_unclassified_percentage = Aquatic_algae_na_count / sumcount_all
    ) %>%
    select(ka, starts_with("Aquatic_algae"))
  
  # --- 2.4 Combine all results for the core ---
  eukaryota_results <- map(c("Fungi", "Metazoa", "Viridiplantae"), calculate_eukaryota_percentages)
  prokaryote_results <- map(c("Bacteria", "Archaea"), calculate_prokaryote_percentages)
  
  all_results_list <- c(
    list(base_counts),
    eukaryota_results, 
    prokaryote_results,
    list(viruses_results), 
    list(aquatic_algae_results)
  ) %>%
    compact() 
  
  final_core_df <- all_results_list %>%
    reduce(full_join, by = "ka") %>%
    rename(Age = ka) %>%
    mutate(Core = core_name, Number = config$id)
  
  if(core_name == "Salmon") {
    final_core_df$Age <- SALMON_NEW_AGES
  }
  
  if(!is.null(config$filter_age)) {
    final_core_df <- final_core_df %>% filter(!Age %in% config$filter_age)
  }
  
  return(final_core_df)
}

#----------------------------------------------------------------#
# 3. EXECUTION: Run the pipeline for all cores
#----------------------------------------------------------------#
all_cores_data <- map2_dfr(names(CORES_CONFIG), CORES_CONFIG, process_core)

write_csv(all_cores_data, file.path(OUTPUT_BASE_PATH, "combine_habitat_percentage_all_cores_refactored.csv"))
message("Successfully saved combined percentage data.")

#----------------------------------------------------------------#
# 4. FINAL CALCULATION: DNA Weight Calculation
#----------------------------------------------------------------#
message("Performing final DNA weight calculations...")

final_processed_data <- all_cores_data %>%
  inner_join(dna_information, by = c("Number", "Age", "Core")) %>%
  mutate(across(c(
    `Total concentration`, `Molecular weight nucleotide`, `Dilution concentration before sequencing 2`,
    `Average read length of sample`, `Library concentration`, `Concentration before Genejet`,
    `DNA`, `Water content`, `Sediment Rate`, `Dry Bulk Density`, `TOC`
  ), as.numeric)) %>%
  mutate(
    Total_concentration_ng_uL = `Total concentration` * 1000,
    Final_DNA_in_Pool = (`Molecular weight nucleotide` * `Dilution concentration before sequencing 2` * `Average read length of sample`) * 0.001,
    Cumulated_Dilution_Factor = Total_concentration_ng_uL / `Dilution concentration before sequencing 2`,
    Cumulated_Concentration_Factor = `Library concentration` / `Concentration before Genejet`,
    Starting_DNA_ng = (Final_DNA_in_Pool * Cumulated_Dilution_Factor) / Cumulated_Concentration_Factor,
    Starting_DNA_g = Starting_DNA_ng * 1e-9,
    Dry_weight = `DNA` * (1 - `Water content` * 0.01),
    Starting_DNA_wt = Starting_DNA_g / Dry_weight,
    ARBULK = `Sediment Rate` * `Dry Bulk Density`,
    BRTOC = (ARBULK * TOC * 0.01)
  )

percentage_columns <- names(all_cores_data)[grepl("_percentage$", names(all_cores_data))]

for (col in percentage_columns) {
  new_col_name <- sub("_percentage$", "_Start_wt", col)
  final_processed_data <- final_processed_data %>%
    mutate(!!new_col_name := Starting_DNA_wt * !!sym(col))
}

write.csv(final_processed_data, file.path(OUTPUT_BASE_PATH, "Processed_DNA_habitat_all_cores_refactored.csv"), row.names = FALSE)

message("Pipeline finished successfully! Final processed data saved.")