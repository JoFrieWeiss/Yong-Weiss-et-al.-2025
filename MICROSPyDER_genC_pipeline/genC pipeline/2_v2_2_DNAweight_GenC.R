# ==============================================================================
# SECOND PART GENC PIPELINE
# ==============================================================================
# please change path for windows in line 24, 25, 27, 35, 38, 395 and 601
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ==============================================================================
# >>> USER CONFIGURATION <<<
# ==============================================================================
# --- B. BACTERIA & ARCHAEA DECISION THRESHOLDS ---
THRESHOLD_AQUATIC     <- -0.5   # Standard:  <= -0.8 is aquatic
THRESHOLD_TERRESTRIAL <- 0.5    # Standard:  >= 0.8 is terrestrial

# ==============================================================================
# 1. LOAD MASTER REFERENCE LISTS
# ==============================================================================
# --- A. PATHS ---
base_dir       <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/"
output_dir <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/output/"
bac_arc_path   <- paste0(base_dir, "code/") 
dna_info_path  <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/Complex_DNA_information.csv"

# --- B. LOAD MASTER LISTS ---
message(">>> Loading Master Lists...")

master_bac <- read_csv(paste0(bac_arc_path, "Master_Bacteria_TAS_Scores.csv"), show_col_types = FALSE)
master_arc <- read_csv(paste0(bac_arc_path, "Master_Archaea_TAS_Scores.csv"), show_col_types = FALSE)

ref_eukaryota <- read_excel("D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/Eukaryota_habitat_new.xlsx") %>%
  rename(family = family) %>% distinct(family, .keep_all = TRUE)

ref_plants_algae <- read_excel("D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/plants_algae_habitat_neu.xlsx") %>%
  rename(family = family) %>% mutate(family = str_trim(family)) %>% distinct(family, .keep_all = TRUE)

# ==============================================================================
# 2. MASTER FUNCTION (With Split Unclassified)
# ==============================================================================

process_core_percentages <- function(core_name, core_data_full) {
  
  message(paste("-----------------------------------------"))
  message(paste("Processing Core:", core_name))
  
  # A. PREPARE CORE TOTALS
  # -----------------------------------------------------
  core_clean_all <- core_data_full %>% filter(ka > 0)
  
  # Gesamtanzahl Reads pro Alter (Nenner für alle Prozente)
  stats_total <- core_clean_all %>%
    group_by(ka) %>%
    summarise(sumcount_all = sum(taxonReads, na.rm = TRUE), .groups="drop")
  
  # Unklassifizierte Reads auf höchster Ebene (Superkingdom NA)
  stats_super_na <- core_clean_all %>%
    filter(is.na(superkingdom)) %>%
    group_by(ka) %>%
    summarise(super_na_sumcount = sum(taxonReads, na.rm = TRUE), .groups="drop")
  
  # Basis-Statistiken zusammenführen
  base_stats <- merge(stats_total, stats_super_na, by="ka", all.x=TRUE) %>%
    mutate(super_na_sumcount = replace_na(super_na_sumcount, 0),
           sumcount_exna = sumcount_all - super_na_sumcount)
  
  
  # B. PROKARYOTES (BACTERIA & ARCHAEA)
  # -----------------------------------------------------
  # Hier bleiben wir bei SPECIES-Level für das Habitat-Matching (TAS Scores)
  process_prokaryote <- function(domain, master_list) {
    
    # 1. Domain Grand Total (ALL Ranks) -> Basis for Taxonomic Unclassified
    dom_total_ALL <- core_clean_all %>% 
      filter(superkingdom == domain) %>% 
      group_by(ka) %>% 
      summarise(dom_sum_reads_ALL = sum(taxonReads), .groups="drop")
    
    # 2. Sum ONLY SPECIES Level -> Basis for Habitat Unclassified
    dom_total_SPECIES <- core_clean_all %>% 
      filter(superkingdom == domain & rank == "species") %>% 
      group_by(ka) %>% 
      summarise(dom_sum_reads_SPECIES = sum(taxonReads), .groups="drop")
    
    # 3. Match with Master List (TAS Scores)
    subset_species <- core_clean_all %>% filter(superkingdom == domain & rank == "species")
    merged <- inner_join(subset_species, master_list, by="species")
    
    # 4. Sum by Habitat
    merged <- merged %>%
      mutate(simple_habitat = case_when(
        tas_score <= THRESHOLD_AQUATIC     ~ "aquatic_1",
        tas_score >= THRESHOLD_TERRESTRIAL ~ "terrestrial_1",
        TRUE ~ "habitat_unclassified_internal"
      ))
    
    
    mid_split <- merged %>%
      filter(simple_habitat == "habitat_unclassified_internal") %>%
      mutate(
        P_aquatic = (THRESHOLD_TERRESTRIAL - tas_score) /
          (THRESHOLD_TERRESTRIAL - THRESHOLD_AQUATIC),
        P_terrestrial = 1 - P_aquatic,
        aquatic_reads = taxonReads * P_aquatic,
        terrestrial_reads = taxonReads * P_terrestrial
      ) %>%
      select(ka, aquatic_reads, terrestrial_reads)
    

    fixed_split <- merged %>%
      filter(simple_habitat %in% c("aquatic_1", "terrestrial_1")) %>%
      mutate(
        aquatic_reads = if_else(simple_habitat == "aquatic_1", taxonReads, 0),
        terrestrial_reads = if_else(simple_habitat == "terrestrial_1", taxonReads, 0)
      ) %>%
      select(ka, aquatic_reads, terrestrial_reads)
    

    grouped <- bind_rows(mid_split, fixed_split) %>%
      group_by(ka) %>%
      summarise(
        aquatic_2 = sum(aquatic_reads),
        terrestrial_2 = sum(terrestrial_reads),
        .groups = "drop"
      )
    
    # Merge everything
    final <- base_stats %>%
      left_join(dom_total_ALL, by="ka") %>%
      left_join(dom_total_SPECIES, by="ka") %>%
      left_join(grouped, by="ka") %>%
      replace(is.na(.), 0)
    
    final %>%
      mutate(
        # A. Taxonomic Unclassified
        # = (All Bacteria - Species Level Bacteria) + Proportion of SuperNA
        # These are reads that weren't taxonomically specific enough (only Genus, Family etc.)
        na_correction = (dom_sum_reads_ALL / sumcount_exna) * super_na_sumcount,
        reads_unc_taxonomic = (dom_sum_reads_ALL - dom_sum_reads_SPECIES) + na_correction,
        
        # B. Habitat Unclassified
        # = (Species Level Bacteria - (Aquatic + Terrestrial))
        # These are reads that ARE a species, but have no habitat label (or are generalists)
        reads_unc_habitat = dom_sum_reads_SPECIES - (aquatic_2 + terrestrial_2),
        
        unclassified_all=+(dom_sum_reads_ALL - (aquatic_2 + terrestrial_2))+na_correction,
        aquatic= aquatic_2+(aquatic_2/(aquatic_2+terrestrial_2))*unclassified_all,
        terrestrial = terrestrial_2+(terrestrial_2/(aquatic_2+terrestrial_2))*unclassified_all,                 
        
        reads_unc_taxonomic_pct=reads_unc_taxonomic/sumcount_all,
        reads_unc_habitat_pct=reads_unc_habitat/sumcount_all,
        # Percentages
        pct_aquatic = aquatic / sumcount_all,
        pct_terrestrial = terrestrial / sumcount_all,
        #pct_unc_tax = reads_unc_taxonomic / sumcount_all,
        #pct_unc_hab = reads_unc_habitat / sumcount_all
        
        aquatic_2_read_pct= aquatic_2,
        terrestrial_2_read_pct= terrestrial_2,
        reads_unc_taxonomic_reads_pct=reads_unc_taxonomic,
        reads_unc_habitat_reads_pct=reads_unc_habitat
        
      ) %>%
      select(ka, pct_aquatic, pct_terrestrial,reads_unc_taxonomic_pct, reads_unc_habitat_pct,na_count = na_correction) %>%
      rename_with(~paste0(domain, "_", .), -ka)
  }
  
  res_bac <- process_prokaryote("Bacteria", master_bac)
  res_arc <- process_prokaryote("Archaea", master_arc)

  
  # C. FUNGI & METAZOA (FAMILY, GENUS & SPECIES)
  # -----------------------------------------------------
  process_general_euk <- function(k_name, ref_db) {
    
    k_total_ALL <- core_clean_all %>% 
      filter(kingdom == k_name) %>% 
      group_by(ka) %>% summarise(k_sum_reads_ALL = sum(taxonReads), .groups="drop")
    
    # Nutzt alles, was mindestens bis zur Familie bestimmt wurde
    subset_family_deeper <- core_clean_all %>% 
      filter(kingdom == k_name & !is.na(family) & family != "")
    
    k_total_FAMILY_DEEPER <- subset_family_deeper %>% 
      group_by(ka) %>% summarise(k_sum_reads_FAMILY = sum(taxonReads), .groups="drop")
    
    merged <- left_join(subset_family_deeper, ref_db, by="family") %>%
      mutate(hab_grp = case_when(
        tolower(habitat) %in% c("freshwater", "marine", "marine/freshwater", "aquatic") ~ "aquatic_1",
        tolower(habitat) %in% c("terrestrial", "soil") ~ "terrestrial_1",
        TRUE ~ "unclassified"
      )) %>%
      group_by(ka, hab_grp) %>% summarise(cnt = sum(taxonReads), .groups="drop") %>%
      pivot_wider(names_from = hab_grp, values_from = cnt, values_fill = 0)
    
    if(!"aquatic_1" %in% names(merged)) merged$aquatic <- 0
    if(!"terrestrial_1" %in% names(merged)) merged$terrestrial <- 0
    
    final <- base_stats %>%
      left_join(k_total_ALL, by="ka") %>%
      left_join(k_total_FAMILY_DEEPER, by="ka") %>%
      left_join(merged, by="ka") %>% replace(is.na(.), 0)
    
    final %>% mutate(
      na_share = (k_sum_reads_ALL / sumcount_exna) * super_na_sumcount,
      reads_unc_tax = (k_sum_reads_ALL - k_sum_reads_FAMILY) + na_share,
      reads_unc_hab = k_sum_reads_FAMILY - (aquatic_1 + terrestrial_1),
      
      unclassified_all=(k_sum_reads_ALL-(aquatic_1 + terrestrial_1))+ na_share,
      
      aquatic=aquatic_1+(aquatic_1/(aquatic_1 + terrestrial_1))*unclassified_all,
      terrestrial=terrestrial_1+(terrestrial_1/(aquatic_1 + terrestrial_1))*unclassified_all,
      
      pct_aqu = aquatic/sumcount_all,
      pct_ter = terrestrial/sumcount_all,
      #pct_utax = reads_unc_tax/sumcount_all,
      #pct_uhab = reads_unc_hab/sumcount_all
    ) %>%
      select(ka, pct_aqu, pct_ter, na_share) %>%
      rename_with(~paste0(k_name, "_", sub("pct_", "", sub("na_share", "na_count", .))), -ka) %>%
      rename_with(~sub("aqu", "aquatic_percentage", .)) %>%
      rename_with(~sub("ter", "terrestrial_percentage", .))

  }
  
  res_fungi <- process_general_euk("Fungi", ref_eukaryota)
  res_metazoa <- process_general_euk("Metazoa", ref_eukaryota)
  
  
  # D. VIRIDIPLANTAE (FAMILY, GENUS & SPECIES)
  # -----------------------------------------------------
  process_plants <- function() {
    
    k_total_ALL <- core_clean_all %>% 
      filter(kingdom == "Viridiplantae") %>% 
      group_by(ka) %>% summarise(k_sum_reads_ALL = sum(taxonReads), .groups="drop")
    
    # Nutzt alles ab Family abwärts
    subset_family_deeper <- core_clean_all %>% 
      filter(kingdom == "Viridiplantae" & !is.na(family) & family != "")
    
    k_total_FAMILY_DEEPER <- subset_family_deeper %>% 
      group_by(ka) %>% summarise(k_sum_reads_FAMILY = sum(taxonReads), .groups="drop")
    
    merged <- left_join(subset_family_deeper, ref_plants_algae, by="family") %>%
      mutate(
        type_cat = case_when(
          tolower(habitat) %in% c("aquatic", "freshwater") ~ "aquatic_plant_1",
          tolower(woody_nonwoody) == "woody" ~ "woody_1",
          tolower(woody_nonwoody) == "non_woody" ~ "non_woody_1",
          TRUE ~ "unclassified"
        )) %>%
      group_by(ka, type_cat) %>% summarise(cnt = sum(taxonReads), .groups="drop") %>%
      pivot_wider(names_from = type_cat, values_from = cnt, values_fill = 0)
    
    for(c in c("aquatic_plant_1","woody_1","non_woody_1")) if(!c %in% names(merged)) merged[[c]]<-0
    
    final <- base_stats %>%
      left_join(k_total_ALL, by="ka") %>%
      left_join(k_total_FAMILY_DEEPER, by="ka") %>%
      left_join(merged, by="ka") %>% replace(is.na(.), 0)
    
    final %>% mutate(
      sum_classified_hab = aquatic_plant_1 + woody_1 + non_woody_1,
      na_share = (k_sum_reads_ALL / sumcount_exna) * super_na_sumcount,
      reads_unc_tax = (k_sum_reads_ALL - k_sum_reads_FAMILY) + na_share,
      reads_unc_hab = k_sum_reads_FAMILY - sum_classified_hab,
      
      unclassified_all=(k_sum_reads_ALL-sum_classified_hab)+ na_share,
      
      aquatic_plant=aquatic_plant_1+((aquatic_plant_1/sum_classified_hab)*unclassified_all),
      woody=woody_1+((woody_1/sum_classified_hab)*unclassified_all),
      non_woody=non_woody_1+((non_woody_1/sum_classified_hab)*unclassified_all),
      
      Viridiplantae_aquatic_percentage=aquatic_plant/sumcount_all,
      Viridiplantae_woody_percentage=woody/sumcount_all,
      Viridiplantae_non_woody_percentage=non_woody/sumcount_all,

      #Viridiplantae_unclassified_taxonomic_percentage = reads_unc_tax/sumcount_all,
      #Viridiplantae_unclassified_habitat_percentage = reads_unc_hab/sumcount_all
    ) %>%
      select(ka, contains("percentage"), Viridiplantae_na_count = na_share)
  }
  res_plants <- process_plants()
  
  
  # E. VIRUSES (SPECIES LEVEL)
  # -----------------------------------------------------
  process_viruses <- function() {
    k_total_ALL <- core_clean_all %>% 
      filter(superkingdom == "Viruses") %>% 
      group_by(ka) %>% summarise(k_sum_reads_ALL = sum(taxonReads), .groups="drop")
    
    k_total_SPECIES <- core_clean_all %>% 
      filter(superkingdom == "Viruses" & rank == "species") %>% 
      group_by(ka) %>% summarise(k_sum_reads_SPECIES = sum(taxonReads), .groups="drop")
    
    final <- base_stats %>% 
      left_join(k_total_ALL, by="ka") %>% 
      left_join(k_total_SPECIES, by="ka") %>%
      replace(is.na(.), 0)
    
    final %>% mutate(
      na_share = (k_sum_reads_ALL / sumcount_exna) * super_na_sumcount,
      reads_unc_tax = (k_sum_reads_ALL - k_sum_reads_SPECIES) + na_share,
      reads_unc_hab = k_sum_reads_SPECIES, 
      
      Viruses_percentage = (k_sum_reads_ALL + na_share)/ sumcount_all,
      #Viruses_unclassified_taxonomic_percentage = reads_unc_tax / sumcount_all,
      #Viruses_unclassified_habitat_percentage = reads_unc_hab / sumcount_all
    ) %>% select(ka, contains("percentage"), Viruses_na_count = na_share)
  }
  res_viruses <- process_viruses()
  
  
  # F. AQUATIC ALGAE (EUKARYOTA NA)
  # -----------------------------------------------------
  process_algae <- function() {
    subset <- core_clean_all %>% filter(is.na(kingdom) & superkingdom == "Eukaryota") 
    total_k <- subset %>% group_by(ka) %>% summarise(k_sum_reads_ALL = sum(taxonReads), .groups="drop")
    
    final <- base_stats %>% left_join(total_k, by="ka") %>% replace(is.na(.), 0)
    
    final %>% mutate(
      na_share = (k_sum_reads_ALL / sumcount_exna) * super_na_sumcount,
      
      Aquatic_algae_percentage = (k_sum_reads_ALL+ na_share)/ sumcount_all,
      
      #Aquatic_algae_unclassified_taxonomic_percentage = na_share / sumcount_all,
      #Aquatic_algae_unclassified_habitat_percentage = 0 
    ) %>% select(ka, contains("percentage"), Aquatic_algae_na_count = na_share)
  }
  res_algae <- process_algae()
  
  
  # G. COMBINE & SAVE
  # -----------------------------------------------------
  list_dfs <- list(res_bac, res_arc, res_fungi, res_plants, res_metazoa, res_viruses, res_algae)
  
  master_table <- Reduce(function(x,y) full_join(x,y, by="ka"), list_dfs) %>%
    arrange(ka) %>%
    rename(Age = ka) 
  
  # Clean column names
  colnames(master_table) <- gsub("_unc_tax", "_unclassified_taxonomic_percentage", colnames(master_table)) %>%
    gsub("_unc_hab", "_unclassified_habitat_percentage", .) %>%
    gsub("_pct_", "_", .) 
  
  final_path <- file.path(output_dir, paste0(core_name, "_FINAL_merged_percentages.csv"))
  write.csv(master_table, final_path, row.names = FALSE)
  
  return(master_table)
}

# ==============================================================================
# 3. EXECUTION (WITH DIAGNOSTICS)
# ==============================================================================

# Ensure the directory exists before starting
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(" Created output directory: ", output_dir)
}

cores_to_process <- list(
  "Ilirney" = "Ilirney",
  "Salmon" = "Salmon",
  "Lele" = "Lele",
  "Lama" = "Lama",
  "Btoko" = "Btoko",
  "Ulu" = "Ulu"
)

message(">>> Starting Processing Loop...")

for (id in names(cores_to_process)) {
  core_obj_name <- cores_to_process[[id]]
  
  if(exists(core_obj_name)) {
    message(paste("Processing core:", core_obj_name))
    
    # Run the function and capture the output
    result <- process_core_percentages(core_obj_name, get(core_obj_name))
    
    # Verify if file was written
    expected_file <- file.path(output_dir, paste0(core_obj_name, "_FINAL_merged_percentages.csv"))
    if(file.exists(expected_file)) {
      message(paste("SUCCESS: File written to", expected_file))
    } else {
      warning(paste("ERROR: Function finished, but NO FILE found at", expected_file))
    }
    
  } else {
    warning(paste("SKIP: Object '", core_obj_name, "' not found in Environment. Did you load the data?"))
  }
}

# ==============================================================================
# SECOND PART: DNA WEIGHT PER TAXON
# ==============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(readxl) 

# ==============================================================================
# 1. DATA PREPARATION (Merge Cores)
# ==============================================================================

# Ensure the path is correct and ends correctly
base_dir_input <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/output"

core_names <- c("Ilirney", "Salmon", "Lele", "Lama", "Btoko", "Ulu")
list_of_dfs <- list()

for (core in core_names) {
  # file.path is safer than paste0 because it handles slashes automatically
  f_path <- file.path(base_dir_input, paste0(core, "_FINAL_merged_percentages.csv"))
  
  # Check for file existence
  if(file.exists(f_path)) {
    message(paste("Found file for core:", core))
    tmp_df <- read_csv(f_path, show_col_types = FALSE) %>%
      mutate(Core = core) 
    list_of_dfs[[core]] <- tmp_df
  } else {
    # This warning will show you exactly where R is looking
    warning(paste("File NOT FOUND at:", f_path))
  }
}

# Combine only if we found any files
if(length(list_of_dfs) > 0) {
  merged_percentages <- bind_rows(list_of_dfs)
  message(paste("Successfully merged", length(list_of_dfs), "core files."))
} else {
  stop("FATAL ERROR: No core files were found! Please check if Part 1 ran correctly and files exist in the output folder.")
}
# B. Load DNA Info
# -------------------------------------------------
DNA_information_unfinal <- read_csv2(dna_info_path, show_col_types = FALSE)


DNA_information <- DNA_information_unfinal %>%
  mutate(
    # Create clean 'Core_Match' column for joining
    Core_Match = case_when(
      str_detect(Core, "Ilirney") ~ "Ilirney",
      str_detect(Core, "Salmon")  ~ "Salmon",
      str_detect(Core, "Lele")    ~ "Lele",
      str_detect(Core, "Lama")    ~ "Lama",
      str_detect(Core, "Btoko")   ~ "Btoko",
      str_detect(Core, "Ulu")     ~ "Ulu",
      TRUE ~ Core # Fallback
    )
  )

# ==============================================================================
# 2. MERGE & CLEANING (FIXED)
# ==============================================================================

# A. Prepare merged_percentages_NEW
merged_percentages_NEW <- merged_percentages %>%
  filter(!(Core == "Ilirney" & round(Age, 4) == 4.400),
         !(Core == "Salmon" & round(Age, 4) == 10.161),
         !(Core == "Lele" & round(Age, 4) == 4.700),
         !(Core == "Ulu" & round(Age, 4) == 35.810)) %>%
  mutate(
    join_Core = trimws(as.character(Core)),
    join_Age  = round(as.numeric(Age), 4)
  )

# B. Prepare DNA_information_NEW (Using your Core_Match logic for joining!)
DNA_information_NEW <- DNA_information %>%
  mutate(
    # Use your already cleaned Core_Match for the join!
    join_Core = trimws(as.character(Core_Match)), 
    join_Age  = round(as.numeric(Age), 4)
  )

# C. THE JOIN
combine_data <- left_join(
  DNA_information_NEW, 
  merged_percentages_NEW, 
  by = c("join_Core", "join_Age"),
  suffix = c(".dna", ".perc") # This prevents name collisions
)

# --- CRITICAL CHECK ---
matched_rows <- sum(!is.na(combine_data$Bacteria_aquatic))
message(paste(">>> DIAGNOSTICS: Successfully matched", matched_rows, "rows with percentage data."))

if(matched_rows == 0) {
  # Let's see why it failed
  print("Example Cores DNA Table:")
  print(unique(DNA_information_NEW$join_Core))
  print("Example Cores Percentage Table:")
  print(unique(merged_percentages_NEW$join_Core))
  stop("FATAL: No matches found between DNA information and Percentages. Check Core names above!")
}

# Cleanup Core/Age columns after join to keep it tidy
combine_data <- combine_data %>%
  mutate(
    Core = coalesce(join_Core),
    Age  = join_Age
  )
# ==============================================================================
# 3. CALCULATIONS
# ==============================================================================

# HERE ARE THE NEW NAMES (adapted to the Split script)
percentage_columns <- c(
  # Bacteria
  "Bacteria_aquatic", 
  "Bacteria_terrestrial", 

  
  # Archaea
  "Archaea_aquatic", 
  "Archaea_terrestrial", 

  
  # Fungi
  "Fungi_aquatic_percentage", 
  "Fungi_terrestrial_percentage", 

  
  # Metazoa
  "Metazoa_aquatic_percentage", 
  "Metazoa_terrestrial_percentage", 

  
  # Viridiplantae (has Woody/Non-Woody)
  "Viridiplantae_aquatic_percentage", 
  "Viridiplantae_woody_percentage", 
  "Viridiplantae_non_woody_percentage", 

  
  # Viruses (has Habitat-Unclassified instead of Ter/Aq)
  "Viruses_percentage", 

  
  # Aquatic Algae (has only Taxonomic Unc.)
  "Aquatic_algae_percentage"

)

# Numeric Conversion
# Include all Percentage-Columns + physical parameters
cols_to_numeric <- c(
  "Total concentration", "Molecular weight nucleotide", "Dilution concentration before sequencing 2",
  "Average read length of sample", "Library concentration", "Concentration before Genejet",
  "DNA", "Water content", "Sediment Rate", "Dry Bulk Density", "TOC",
  percentage_columns 
)

# Convert only existing columns
existing_cols <- intersect(cols_to_numeric, names(combine_data))

combine_data <- combine_data %>%
  mutate(across(all_of(existing_cols), ~as.numeric(as.character(.))))

# --- Basis Calculations (DNA Weight) ---
combine_data_calc <- combine_data %>%
  mutate(
    # 1. Base Values
    Total_concentration_calc = `Total concentration` * 1000,
    
    Final_DNA_in_Pool = (`Molecular weight nucleotide` * `Dilution concentration before sequencing 2` * `Average read length of sample`) * 0.001,
    
    Cumulated_Dilution_Factor = Total_concentration_calc / `Dilution concentration before sequencing 2`,
    
    Cumulated_Concentration_Factor = `Library concentration` / `Concentration before Genejet`,
    
    # 2. Starting DNA
    Starting_DNA_ng = (Final_DNA_in_Pool * Cumulated_Dilution_Factor) / Cumulated_Concentration_Factor,
    Starting_DNA_g = Starting_DNA_ng * 1e-9,
    
    # 3. Dry Weight Normalization
    Dry_weight = `DNA` * (1 - `Water content` * 0.01),
    
    # THE MOST IMPORTANT: g DNA per g Sediment (Total)
    Starting_DNA_wt = Starting_DNA_g / Dry_weight,
    
    # Sediment Accumulation
    ARBULK = `Sediment Rate` * `Dry Bulk Density`,
    BRTOC = (ARBULK * TOC * 0.01)
  )

# --- Loop for Specific Weights ---
# Calculates e.g. "Bacteria_aquatic_percentage_Start_wt"
for (col in percentage_columns) {
  new_col_name <- paste0(col, "_Start_wt")
  
  if (col %in% names(combine_data_calc)) {
    # Calculation: Proportion * Total Weight
    combine_data_calc[[new_col_name]] <- combine_data_calc[[col]] * combine_data_calc$Starting_DNA_wt
  } else {
    # Optional warning; if Algae lacks a column, make it NA
    # (With the new lists everything should match)
    # warning(paste("Column missing for weight calc:", col)) 
    combine_data_calc[[new_col_name]] <- NA
  }
}

# ==============================================================================
# 4. SAVE
# ==============================================================================

# Adjust path if you want to save a new version
write_csv(combine_data_calc, "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/output/Processed_DNA_habitat_percentage_SPLIT_UNCLASSIFIED.csv")

message("Calculation complete! File with split unclassified values saved.")
