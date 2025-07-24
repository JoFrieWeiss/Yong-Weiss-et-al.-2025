# =============================================================================
#
# COMPLETE UNIFIED ANALYSIS PIPELINE FOR SEDADNA HABITAT ANALYSIS
#
# DESCRIPTION:
# This single, master script performs the entire analysis workflow.
# It is designed to be run from top to bottom.
#
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 1: CONFIGURATION (All user settings are here)
# -----------------------------------------------------------------------------

# --- 1.1: Load Libraries ---
library(tidyverse)
library(readr)
library(readxl)
library(stringr)
library(purrr)

# --- 1.2: Define Main Paths ---
BASE_PATH <- "/Volumes/projects/biodiv/user/ziyong"

# Directories for outputs
INTERMEDIATE_DIR <- file.path(BASE_PATH, "A_aquatic_terrestrial/1_Prokaryote_Classification_Output")
FINAL_OUTPUT_DIR   <- file.path(BASE_PATH, "A_aquatic_terrestrial/2_Final_Analysis_Output")

# --- 1.3: Configure Input File Paths ---
# Using the paths you provided
core_rdata_paths <- c(
  Btoko   = file.path(BASE_PATH, "A_shotgun_data/btoko_lca.v1_MergedData.RData"),
  Ilirney = file.path(BASE_PATH, "A_shotgun_data/ilirney_lca.v1_MergedData.RData"),
  Lama    = file.path(BASE_PATH, "A_shotgun_data/Lama_lca_MergedData.RData"),
  Lele    = file.path(BASE_PATH, "A_shotgun_data/lele_lca.v1_MergedData.RData"),
  Salmon  = file.path(BASE_PATH, "A_shotgun_data/salmon_lca.v2_MergedData.RData"),
  Ulu     = file.path(BASE_PATH, "A_shotgun_data/ulu_lca.v1_MergedData.RData")
)

bacteria_habitat_paths <- c(
  Ilirney = file.path(BASE_PATH, "1_ilirney Lake/species_list/updated_Bacteria_Ilirney_list.csv"),
  Salmon  = file.path(BASE_PATH, "2_Salmon Lake/species_list/updated_Bacteria_Salmon_list.csv"),
  Lele    = file.path(BASE_PATH, "3_Lele Lake/species_list/updated_Bacteria_Lele_list.csv"),
  Btoko   = file.path(BASE_PATH, "5_Btoko Lake/species_list/updated_bacteria_Btoko_list.csv"),
  Lama    = file.path(BASE_PATH, "4_Lama Lake/species_list/updated_Bacteria_Lama_list.csv"),
  Ulu     = file.path(BASE_PATH, "6_ulu Lake/species_list/updated_Bacteria_ulu_list.csv")
)

archaea_habitat_paths <- c(
  Ilirney = file.path(BASE_PATH, "1_ilirney Lake/species_list/updated_Archaea_Ilirney_list.csv"),
  Salmon  = file.path(BASE_PATH, "2_Salmon Lake/species_list/updated_Archaea_Salmon_list.csv"),
  Lele    = file.path(BASE_PATH, "3_Lele Lake/species_list/updated_Archaea_Lele_list.csv"),
  Btoko   = file.path(BASE_PATH, "5_Btoko Lake/species_list/updated_archaea_Btoko_list.csv"),
  Lama    = file.path(BASE_PATH, "4_Lama Lake/species_list/updated_Archaea_Lama_list.csv"),
  Ulu     = file.path(BASE_PATH, "6_ulu Lake/species_list/updated_Archaea_ulu_list.csv")
)

EUKARYOTA_HABITAT_DB_PATH <- "/Users/josefineweiss/Downloads/final_results_eukaryota_excel-4.xlsx"
PLANT_FAMILY_DB_PATH <- file.path(BASE_PATH, "A_Eukaryota_Habitat/Plant_Family_Type.csv")


# -----------------------------------------------------------------------------
# SECTION 2: HELPER FUNCTIONS (The complete "Toolbox" for the pipeline)
# -----------------------------------------------------------------------------

preprocess_core_data <- function(core_data) {
  data_filtered <- core_data %>% filter(ka > 0)
  total_reads_df <- data_filtered %>% group_by(ka) %>% summarise(total_reads = sum(taxonReads, na.rm = TRUE), .groups = "drop")
  superkingdom_na_reads_df <- data_filtered %>% filter(is.na(superkingdom)) %>% group_by(ka) %>% summarise(super_na_reads = sum(taxonReads, na.rm = TRUE), .groups = "drop")
  base_metrics <- total_reads_df %>%
    left_join(superkingdom_na_reads_df, by = "ka") %>%
    mutate(super_na_reads = replace_na(super_na_reads, 0), reads_without_na = total_reads - super_na_reads)
  return(list(base_metrics = base_metrics, data_filtered = data_filtered))
}

assign_eukaryota_habitats <- function(eukaryota_data, habitat_db) {
  habitat_db_renamed <- habitat_db %>% rename(species = Species)
  combined_data <- left_join(eukaryota_data, habitat_db_renamed, by = "species")
  refined_data <- combined_data %>%
    mutate(new_habitat = case_when(
      phylum == "Streptophyta" & !family %in% c("Butomaceae", "Characeae", "Closteriaceae", "Desmidiaceae", "Elatinaceae", "Equisetaceae", "Hydrocharitaceae", "Juncaceae", "Nymphaeaceae", "Lentibulariaceae", "Potamogetonaceae", "Typhaceae", "Juncaginaceae") ~ "terrestrial",
      order == "Lepidoptera" ~ "terrestrial", family == "Potamogetonaceae" ~ "freshwater", species == "Allium schoenoprasum" ~ "terrestrial",
      species %in% c("Equisetum hyemale", "Luzula luzuloides", "Sparganium erectum", "Triglochin maritima") ~ "freshwater/terrestrial",
      species == "Penaeus chinensis" | order == "Octopoda" ~ "marine", family == "Cordylidae" ~ "terrestrial", genus == "Cherax" ~ "freshwater",
      species %in% c("Boleophthalmus pectinirostris", "Salvelinus fontinalis") ~ "freshwater", TRUE ~ Habitat
    ))
  return(refined_data)
}

habitat_decision_making <- function(data) {
  data %>%
    mutate(habitat = case_when(
      str_detect(isolation_categories, regex("#(freshwater|lake|river|pond|marine)", ignore_case = TRUE)) ~ "aquatic",
      str_detect(isolation_categories, regex("#(terrestrial|soil|sediment|mud|dust)", ignore_case = TRUE)) ~ "terrestrial",
      ratio_terrestrial > 0.2 & ratio_aquatic > 0.2 ~ "aquatic/terrestrial",
      ratio_terrestrial > 0.2 ~ "terrestrial", ratio_aquatic > 0.2 ~ "aquatic", TRUE ~ "unclassified"
    ))
}

prepare_and_classify_prokaryotes <- function(core_name, kingdom_name, rdata_path, habitat_csv_path, output_dir) {
  habitat_info <- read_csv(habitat_csv_path, show_col_types = FALSE)
  temp_env <- new.env(); load(rdata_path, envir = temp_env); merged_data <- temp_env$MergedData
  kingdom_reads <- merged_data %>% filter(superkingdom == !!kingdom_name, rank == "species", ka > 0)
  habitat_info_processed <- habitat_info %>%
    mutate(Total_sample_number = aquatic_samples + soil_samples + animal_samples + plant_samples,
           ratio_terrestrial = soil_samples / Total_sample_number, ratio_aquatic = aquatic_samples / Total_sample_number) %>%
    mutate(across(starts_with("ratio_"), ~replace_na(., 0)))
  combined_data <- merge(kingdom_reads, habitat_info_processed, by = "species")
  classified_data <- habitat_decision_making(combined_data)
  output_filename <- file.path(output_dir, paste0(core_name, "_classified_", kingdom_name, ".csv"))
  final_export <- classified_data %>% select(ka, species, taxonReads, habitat)
  write_csv(final_export, output_filename)
}

calculate_eukaryote_habitat_percentages <- function(kingdom_name, base_metrics, all_filtered_data, annotated_eukaryotes) {
  kingdom_summary <- all_filtered_data %>% filter(kingdom == !!kingdom_name) %>% group_by(ka) %>% summarise(kingdom_total = sum(taxonReads, na.rm = TRUE))
  kingdom_annotated <- annotated_eukaryotes %>% filter(kingdom == !!kingdom_name)
  
  kingdom_grouped <- kingdom_annotated %>%
    group_by(ka, new_habitat) %>%
    summarise(sumcount = sum(taxonReads, na.rm = TRUE), .groups = "drop") %>%
    mutate(habitat_group = case_when(
      new_habitat %in% c("freshwater", "marine", "marine/freshwater") ~ "aquatic",
      new_habitat == "terrestrial" ~ "terrestrial", TRUE ~ "other"
    )) %>%
    group_by(ka, habitat_group) %>%
    summarise(total_in_group = sum(sumcount, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = habitat_group, values_from = total_in_group, values_fill = 0)
  
  # *** KORREKTUR IST HIER ***
  # Die gesamte Berechnungskette wird jetzt der Variable 'result' zugewiesen.
  result <- reduce(list(base_metrics, kingdom_summary, kingdom_grouped), full_join, by = "ka") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    mutate(
      ter_aqu_total = rowSums(across(any_of(c("terrestrial", "aquatic"))), na.rm = TRUE),
      kingdom_na_estimate = (kingdom_total / reads_without_na) * super_na_reads,
      kingdom_unclassified = (kingdom_total - ter_aqu_total) + kingdom_na_estimate,
      perc_aquatic = aquatic / total_reads,
      perc_terrestrial = terrestrial / total_reads,
      perc_unclassified = kingdom_unclassified / total_reads,
      perc_total = kingdom_total / total_reads
    ) %>%
    select(ka, starts_with("perc_"))
  
  names(result) <- paste(kingdom_name, names(result), sep = "_")
  return(result %>% rename(ka = !!paste0(kingdom_name, "_ka")))
}

calculate_simple_percentage <- function(group_name, preprocessed_info) {
  base_metrics <- preprocessed_info$base_metrics
  data_filtered <- preprocessed_info$data_filtered
  
  if (group_name == "Aquatic_algae") {
    group_data <- data_filtered %>% filter(superkingdom == "Eukaryota", is.na(kingdom))
  } else {
    group_data <- data_filtered %>% filter(superkingdom == !!group_name)
  }
  
  total_group_reads <- group_data %>% group_by(ka) %>% summarise(sumcount_group_all = sum(taxonReads, na.rm = TRUE), .groups = "drop")
  
  final_df <- base_metrics %>%
    left_join(total_group_reads, by = "ka") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    mutate(
      na_count = (sumcount_group_all / reads_without_na) * super_na_reads,
      perc_total = sumcount_group_all / total_reads,
      perc_unclassified = na_count / total_reads
    ) %>%
    select(ka, perc_total, perc_unclassified)
  
  names(final_df) <- paste(group_name, names(final_df), sep = "_")
  return(final_df %>% rename(ka = !!paste0(group_name, "_ka")))
}

# =============================================================================
# SECTION 3: MAIN WORKFLOW EXECUTION
# =============================================================================

# --- Create Output Directories ---
if (!dir.exists(INTERMEDIATE_DIR)) dir.create(INTERMEDIATE_DIR, recursive = TRUE)
if (!dir.exists(FINAL_OUTPUT_DIR)) dir.create(FINAL_OUTPUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# WORKFLOW STEP 1: Prepare and Classify Prokaryote Habitats
# -----------------------------------------------------------------------------
cat("### STARTING WORKFLOW STEP 1: Preparing Prokaryote Data... ###\n\n")

for (core in names(core_rdata_paths)) {
  for (kingdom in c("Bacteria", "Archaea")) {
    habitat_path <- if (kingdom == "Bacteria") bacteria_habitat_paths[core] else archaea_habitat_paths[core]
    if (is.na(habitat_path) || !file.exists(habitat_path)) {
      warning(paste("Skipping:", core, "-", kingdom, "-> Input CSV not found or configured."))
      next
    }
    cat(paste(" -> Classifying", kingdom, "for", core, "core...\n"))
    prepare_and_classify_prokaryotes(
      core_name = core, kingdom_name = kingdom,
      rdata_path = core_rdata_paths[core], habitat_csv_path = habitat_path,
      output_dir = INTERMEDIATE_DIR
    )
  }
}
cat("\n### WORKFLOW STEP 1 COMPLETE ###\n\n")

# -----------------------------------------------------------------------------
# WORKFLOW STEP 2: Calculate Final Percentages for All Kingdoms
# -----------------------------------------------------------------------------
cat("### STARTING WORKFLOW STEP 2: Calculating Final Percentages... ###\n\n")

eukaryota_habitat_db <- read_excel(EUKARYOTA_HABITAT_DB_PATH)
plant_family_db <- read_csv(PLANT_FAMILY_DB_PATH, show_col_types = FALSE)
all_core_results <- list()

for (core_name in names(core_rdata_paths)) {
  cat(paste(" -> Processing all kingdoms for:", core_name, "core...\n"))
  load(core_rdata_paths[core_name])
  preprocessed_info <- preprocess_core_data(MergedData)
  current_core_results <- list()
  
  eukaryota_reads <- preprocessed_info$data_filtered %>% filter(superkingdom == "Eukaryota", rank == "species")
  eukaryota_with_habitats <- assign_eukaryota_habitats(eukaryota_reads, eukaryota_habitat_db)
  
  current_core_results$Fungi <- calculate_eukaryote_habitat_percentages("Fungi", preprocessed_info$base_metrics, preprocessed_info$data_filtered, eukaryota_with_habitats)
  current_core_results$Viridiplantae <- calculate_eukaryote_habitat_percentages("Viridiplantae", preprocessed_info$base_metrics, preprocessed_info$data_filtered, eukaryota_with_habitats)
  current_core_results$Metazoa <- calculate_eukaryote_habitat_percentages("Metazoa", preprocessed_info$base_metrics, preprocessed_info$data_filtered, eukaryota_with_habitats)
  current_core_results$Viruses <- calculate_simple_percentage("Viruses", preprocessed_info)
  current_core_results$Aquatic_algae <- calculate_simple_percentage("Aquatic_algae", preprocessed_info)
  
  for (kingdom in c("Bacteria", "Archaea")) {
    prokaryote_file <- file.path(INTERMEDIATE_DIR, paste0(core_name, "_classified_", kingdom, ".csv"))
    if (file.exists(prokaryote_file)) {
      classified_prokaryotes <- read_csv(prokaryote_file, show_col_types = FALSE)
      prokaryote_habitat_counts <- classified_prokaryotes %>% group_by(ka, habitat) %>% summarise(read_count = sum(taxonReads, na.rm = TRUE), .groups = "drop") %>% pivot_wider(names_from = habitat, values_from = read_count, values_fill = 0)
      total_kingdom_reads <- preprocessed_info$data_filtered %>% filter(superkingdom == !!kingdom) %>% group_by(ka) %>% summarise(sumcount_kingdom_all = sum(taxonReads, na.rm = TRUE), .groups = "drop")
      prokaryote_summary <- preprocessed_info$base_metrics %>%
        left_join(total_kingdom_reads, by = "ka") %>% left_join(prokaryote_habitat_counts, by = "ka") %>%
        mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
        mutate(
          na_count = (sumcount_kingdom_all / reads_without_na) * super_na_reads,
          classified_sum = rowSums(across(any_of(c("aquatic", "terrestrial"))), na.rm = TRUE),
          unclassified_count = (sumcount_kingdom_all - classified_sum) + na_count,
          perc_aquatic = aquatic / total_reads,
          perc_terrestrial = terrestrial / total_reads,
          perc_unclassified = unclassified_count / total_reads
        ) %>% select(ka, starts_with("perc_"))
      names(prokaryote_summary) <- paste(kingdom, names(prokaryote_summary), sep = "_")
      current_core_results[[kingdom]] <- prokaryote_summary %>% rename(ka = !!paste0(kingdom, "_ka"))
    }
  }
  
  combined_core_df <- Reduce(function(x, y) full_join(x, y, by = "ka"), current_core_results)
  all_core_results[[core_name]] <- combined_core_df
  rm(MergedData)
  cat("  -> Finished processing", core_name, "\n")
}
cat("\n### WORKFLOW STEP 2 COMPLETE ###\n\n")

# -----------------------------------------------------------------------------
# WORKFLOW STEP 3: Combine All Cores and Save Final Report
# -----------------------------------------------------------------------------

cat("### STARTING WORKFLOW STEP 3: Aggregating Final Report... ###\n")

final_report <- bind_rows(all_core_results, .id = "Core")
final_report <- final_report %>%
  rename(Age = ka) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  select(Core, Age, everything())

final_output_path <- file.path(FINAL_OUTPUT_DIR, "combine_habitat_percentage_aqu_ter_fre_unclassfied_new.csv")
write_csv(final_report, final_output_path)

cat(paste("ANALYSIS PIPELINE FINISHED! Final report saved to:", final_output_path))

