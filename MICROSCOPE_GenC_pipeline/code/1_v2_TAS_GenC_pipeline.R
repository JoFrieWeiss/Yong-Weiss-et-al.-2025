# ==============================================================================
# FIRST PART GENC PIPELINE
# ==============================================================================
# please change path for windows in line 19
# ==============================================================================

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyr)
if(!require(ggrepel)) install.packages("ggrepel")
library(ggrepel)

# ==============================================================================
# LOAD SHOTGUN DATA——new
# ==============================================================================

# Assumes data objects (Ilirney, Btoko, etc.) are already loaded in the environment
# Example placeholders:
Ilirney <- MergedData
Salmon <- MergedData
Lele <- MergedData
Lama <- MergedData
Btoko <- MergedData
Ulu <- MergedData

# ==============================================================================
# 1. SETUP & CLEANUP OF OLD FILES
# ==============================================================================
base_dir <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026"
base_dir_code <- file.path(base_dir,"code/")
path_bac_master <- file.path(base_dir_code, "Master_Bacteria_TAS_Scores.csv")
path_arc_master <- file.path(base_dir_code, "Master_Archaea_TAS_Scores.csv")

if(file.exists(path_bac_master)) {
  file.remove(path_bac_master)
  message("Old Master_Bacteria_TAS_Scores.csv deleted.")
}

# ==============================================================================
# 2. FUNCTION: CLEAN BUILD (WITH STRINGENT ARTIFACT/ZOMBIE FILTER)
# ==============================================================================
create_clean_master <- function(filepath, target_domain) {
  if(is.null(filepath) || !file.exists(filepath)) return(NULL)
  
  message(paste("\n>>> PROCESSING:", target_domain))
  raw_data <- read_csv2(filepath, show_col_types = FALSE)
  
  # Standardize column names
  weird_name <- "unique(Master_Bacteria_Raw$species)"
  if(weird_name %in% names(raw_data)) raw_data <- raw_data %>% rename(species = all_of(weird_name))
  if("Name" %in% names(raw_data)) raw_data <- raw_data %>% rename(species = Name)
  
  # Habitat-specific column lists
  cols_aquatic <- c("Total_Aquatic", "Aquatic_Marine", "Aquatic_Ocean", "Aquatic_Sea", "Aquatic_Sediment", "Aquatic_Estuary", "Aquatic_Waste_Water", "Aquatic_River", "Aquatic_Ice", "Aquatic_Lake", "Aquatic_Unknown", "Aquatic_Groundwater", "Aquatic_Reservoir", "Aquatic_Brine", "Animal_Shark", "Animal_Fish", "Animal_Whale", "Animal_Dolphin", "Animal_Tadpole")
  cols_terrestrial <- c("Total_Soil", "Total_Plant", "Soil_Desert", "Soil_Farm", "Soil_Peatland", "Soil_Field", "Soil_Forest", "Soil_Unknown", "Soil_Tundra", "Soil_Agricultural", "Soil_Paddy", "Soil_Shrub", "Plant_Rhizosphere", "Plant_Leaf", "Plant_Unknown", "Plant_Flower", "Plant_Seed", "Plant_Stem", "Plant_Wood", "Plant_Sprout", "Animal_Human", "Animal_Goat", "Animal_Bird", "Animal_Cattle", "Animal_Mouse", "Animal_Swallow", "Animal_Sheep", "Animal_Dog", "Animal_Rat", "Animal_Horse", "Animal_Pig", "Animal_Termite", "Animal_Fly", "Animal_Bee", "Animal_Mosquito", "Animal_Baboon", "Animal_Cat", "Animal_Tick", "Animal_Bat", "Animal_Macaque", "Animal_Bumblebee", "Animal_Chimpanzee", "Animal_Sparrow", "Animal_Gorilla", "Animal_Roe", "Animal_Pigeon", "Animal_Seagull", "Animal_Hamster", "Animal_Fruitfly")
  
  # Helper to sum relevant columns and convert numeric formats
  sum_cols <- function(df, col_list) {
    valid <- intersect(names(df), col_list)
    if(length(valid)==0) return(rep(0, nrow(df)))
    df %>% select(all_of(valid)) %>%
      mutate(across(everything(), ~ as.numeric(gsub(",", ".", as.character(.))))) %>%
      mutate(across(everything(), ~ replace_na(., 0))) %>%
      rowSums(na.rm=TRUE)
  }
  
  clean_df <- raw_data %>%
    mutate(
      Val_Aqua = sum_cols(., cols_aquatic),
      Val_Terr = sum_cols(., cols_terrestrial),
      Total_Sum = Val_Aqua + Val_Terr
    ) %>%
    filter(Total_Sum > 0) %>%
    mutate(
      tas_score = (Val_Terr - Val_Aqua) / Total_Sum
    ) %>%
    # --- ARTIFACT FILTER (Excluding specific artificial value 0.3461187) ---
    filter(abs(tas_score - 0.3461187) > 0.0000001) %>%
    
    mutate(
      habitat_label = case_when(
        tas_score >= 0.8  ~ "Terrestrial Specialist",
        tas_score > 0.1   ~ "Mostly Terrestrial",
        tas_score <= -0.8 ~ "Aquatic Specialist",
        tas_score < -0.1  ~ "Mostly Aquatic",
        TRUE ~ "Generalist"
      ),
      total_n = 9999, r_terr = ifelse(tas_score > 0, 1, 0), r_aqua = ifelse(tas_score < 0, 1, 0),
      r_anim = 0, r_plant = 0, tas_score_raw = tas_score
    ) %>%
    select(species, tas_score, habitat_label, total_n, r_terr, r_aqua, r_anim, r_plant, tas_score_raw) %>%
    distinct(species, .keep_all = TRUE)
  
  return(clean_df)
}

# ==============================================================================
# 3. CALCULATE DATA & SAVE MASTER LISTS
# ==============================================================================
file_bacteria <- file.path(base_dir_code, "MicrobeAtlas_Full_Results_ALLBACTERIA.csv")
Master_Bacteria_TAS <- create_clean_master(file_bacteria, "Bacteria")

if(!is.null(Master_Bacteria_TAS)) {
  write_csv(Master_Bacteria_TAS, path_bac_master)
  message("Cleaned Bacteria list saved. Artifacts removed.")
}

# Load Archaea (if available)
if(file.exists(path_arc_master)) {
  Master_Archaea_TAS <- read_csv(path_arc_master, show_col_types = FALSE)
} else {
  Master_Archaea_TAS <- tibble(species = character())
}

# ==============================================================================
# 4. FINAL PLOT (HISTOGRAM WITHOUT ARTIFACTS)
# ==============================================================================
message("--- GENERATING PLOT (ARTIFACT-FREE) ---")

core_names_list <- c("Ilirney", "Salmon", "Lele", "Lama", "Btoko", "Ulu")
all_reads_list <- list()

for (core_name in core_names_list) {
  if(exists(core_name)) {
    df <- get(core_name)
    if("taxonReads" %in% names(df)) {
      temp <- df %>% 
        filter(rank == "species") %>%
        group_by(species, superkingdom) %>%
        summarise(Reads = sum(taxonReads), .groups = "drop")
      all_reads_list[[core_name]] <- temp
    }
  }
}

global_reads <- bind_rows(all_reads_list) %>%
  group_by(species, superkingdom) %>%
  summarise(Total_Reads = sum(Reads), .groups = "drop")

reads_bac <- global_reads %>%
  filter(superkingdom == "Bacteria") %>%
  inner_join(Master_Bacteria_TAS, by = "species") %>%
  mutate(Domain_Label = "Bacteria")

reads_arch <- global_reads %>%
  filter(superkingdom == "Archaea") %>%
  inner_join(Master_Archaea_TAS, by = "species") %>%
  mutate(Domain_Label = "Archaea")

plot_data <- bind_rows(reads_bac, reads_arch)

if(nrow(plot_data) > 0) {
  get_top_labels <- function(data, n_top = 1) {
    data %>%
      mutate(Zone = cut(tas_score, breaks = c(-1.1, -0.6, -0.2, 0.2, 0.6, 1.1), labels = c("Aq+","Aq","Gen","Terr","Terr+"))) %>%
      group_by(Domain_Label, Zone) %>%
      arrange(desc(Total_Reads)) %>%
      slice_head(n = n_top) %>% ungroup()
  }
  top_labels <- get_top_labels(plot_data, n_top = 2)
  
  n_bac <- nrow(reads_bac); n_arch <- nrow(reads_arch)
  plot_data <- plot_data %>% mutate(Facet_Title = ifelse(Domain_Label=="Bacteria", paste0("Bacteria (n=",n_bac,")"), paste0("Archaea (n=",n_arch,")")))
  top_labels <- top_labels %>% mutate(Facet_Title = ifelse(Domain_Label=="Bacteria", paste0("Bacteria (n=",n_bac,")"), paste0("Archaea (n=",n_arch,")")))
  
  p_clean <- ggplot(plot_data, aes(x = tas_score)) +
    geom_histogram(aes(fill = ..x..), binwidth = 0.05, center = 0.025, color = "white", size=0.1) +
    geom_label_repel(data = top_labels, aes(y = 10, label = species), size = 2.5, fontface = "italic", max.overlaps = 30, nudge_y = 50) +
    facet_wrap(~Facet_Title, scales = "free_y") +
    scale_fill_gradient2(low = "#0072B2", mid = "#eeeeee", high = "#D55E00", midpoint = 0, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title = "Global TAS (Artifacts Filtered)", x = "TAS Score", y = "Count") +
    theme_minimal()
  
  ggsave(paste0(base_dir_code,"TAS_over_all_cores.png"), p_clean, width = 14, height = 8, bg = "white")
  print(p_clean)
}

# ==============================================================================
# 5. ECO-DRIVER PLOT (ARTIFACT-FREE)
# ==============================================================================
message("--- GENERATING ECO-DRIVER PLOT ---")

# Reload raw data for driver analysis and join with filtered species
habitat_profile <- Master_Bacteria_TAS %>%
  inner_join(
    read_csv2(file_bacteria, show_col_types = FALSE) %>% 
      rename(species = 1) %>% 
      mutate(across(c(Total_Animal, Total_Plant, Total_Soil, Total_Aquatic), ~ as.numeric(gsub(",", ".", .)))),
    by = "species"
  ) %>%
  mutate(TAS_Bin = round(tas_score, 1)) %>%
  group_by(TAS_Bin) %>%
  summarise(
    Avg_Soil = mean(Total_Soil, na.rm = TRUE),
    Avg_Water = mean(Total_Aquatic, na.rm = TRUE),
    Avg_Animal = mean(Total_Animal, na.rm = TRUE),
    Avg_Plant = mean(Total_Plant, na.rm = TRUE),
    n_species = n()
  ) %>%
  pivot_longer(cols = starts_with("Avg"), names_to = "Habitat", values_to = "Density")

p_profile <- ggplot(habitat_profile, aes(x = TAS_Bin, y = Density, color = Habitat, group = Habitat)) +
  geom_line(size = 1.2) +
  geom_point(aes(size = n_species), alpha = 0.6) +
  scale_color_manual(values = c("Avg_Soil"="#2ca02c", "Avg_Water"="#1f77b4", "Avg_Animal"="#d62728", "Avg_Plant"="#ff7f0e")) +
  labs(title = "Ecological Drivers (Artifact-Free)", subtitle = "Showing which habitats dominate at which score level", x = "TAS Score", y = "Average Frequency") +
  theme_minimal()

print(p_profile)

# ==============================================================================
# 5b. ANIMAL-SPECIFIC DRIVER PLOT
# ==============================================================================
message("--- GENERATING ANIMAL-SPECIFIC DRIVER PLOT ---")

# Define detailed animal subgroups (ensure these match CSV headers exactly)
animal_subgroups <- c("Animal_Human", "Animal_Fish", "Animal_Cattle", 
                      "Animal_Dog", "Animal_Mouse", "Animal_Bird")

# Load raw data and extract specific animal host frequencies
animal_profile <- Master_Bacteria_TAS %>%
  inner_join(
    read_csv2(file_bacteria, show_col_types = FALSE) %>% 
      rename(species = 1) %>% 
      # Convert all relevant animal columns to numeric
      mutate(across(all_of(intersect(names(.), c(animal_subgroups, "Total_Animal"))), 
                    ~ as.numeric(gsub(",", ".", .)))),
    by = "species"
  ) %>%
  mutate(TAS_Bin = round(tas_score, 1)) %>%
  group_by(TAS_Bin) %>%
  summarise(
    # Calculate mean frequency for each animal group per bin
    across(all_of(intersect(names(.), animal_subgroups)), mean, na.rm = TRUE),
    n_species = n()
  ) %>%
  # Reshape for ggplot
  pivot_longer(cols = starts_with("Animal_"), names_to = "Animal_Group", values_to = "Density")

# Plot: Specific animal host drivers
p_animal_drivers <- ggplot(animal_profile, aes(x = TAS_Bin, y = Density, color = Animal_Group, group = Animal_Group)) +
  geom_line(size = 1.1, alpha = 0.8) +
  geom_point(aes(size = n_species), alpha = 0.5) +
  scale_color_brewer(palette = "Set1") + # Distinct color palette
  labs(
    title = "Animal-Specific Drivers of TAS Score",
    subtitle = "Which animal hosts dominate at which score level?",
    x = "TAS Score (-1 = Aquatic, +1 = Terrestrial)",
    y = "Average Frequency within Host Group",
    color = "Animal Host"
  ) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

print(p_animal_drivers)

# Save output
ggsave(paste0(base_dir_code,"ANIMAL_SPECIFIC_DRIVERS.png"), 
       p_animal_drivers, width = 12, height = 7, bg = "white")

message("--- PART 1 OF GENC PIPELINE DONE ---")
