# =====================================================
# Title: Multivariate Analyses of Aquatic-Terrestrial Transition Data
# Author: [Zijuan Yong]
# Date: [2025-07-15]
# Description: Performs PCA, RDA, and GLMM on sedimentary DNA biomass and environmental data.
# =====================================================

# =====================
# 1. Load Required Libraries
# =====================
library(vegan)        # PCA & RDA
library(lmerTest)     # GLMM with lmer() and p-values
library(ggeffects)    # Marginal effects
library(dplyr)        # Data manipulation

# =====================
# 2. Set Paths
# =====================
data_dir <- "Data"
figure_dir <- "Figure"

# =====================
# 3. Redundancy Analysis (RDA)
# =====================

## Load and preprocess biomass data
biomass_data <- read.csv(paste0(data_dir, "Calculate_all_information_result_unclassfied.csv"))
biomass_data[is.na(biomass_data)] <- 0

## Assign timeslices
timeslices <- seq(0, 60, 2)
assign_timeslice <- function(Age) {
  for (i in 1:(length(timeslices) - 1)) {
    if (timeslices[i] <= Age & Age < timeslices[i + 1]) {
      return(paste0(timeslices[i], "-", timeslices[i + 1]))
    }
  }
  return(paste0(timeslices[length(timeslices)], "+"))
}
biomass_data$timeslice <- sapply(biomass_data$Age, assign_timeslice)
biomass_data$Startage <- as.integer(sub("-.*", "", biomass_data$timeslice))

## Assign lake group
biomass_data <- biomass_data %>%
  mutate(core_name = case_when(
    grepl("Ilirney", Core, ignore.case = TRUE) ~ "Ilirney",
    grepl("Lele", Core, ignore.case = TRUE) ~ "Lele",
    grepl("Lama", Core, ignore.case = TRUE) ~ "Lama",
    grepl("Btoko", Core, ignore.case = TRUE) ~ "BToko",
    grepl("Ulu", Core, ignore.case = TRUE) ~ "Ulu",
    grepl("Salmon", Core, ignore.case = TRUE) ~ "Salmon",
    TRUE ~ NA_character_
  ),
  Lake_group = factor(core_name))

## Species data (fourth-root transformed)
species_data <- biomass_data %>%
  select(starts_with("median_biomass")) %>%
  rename_with(~ gsub("median_biomass_|_percent", "", .x)) %>%
  rename_with(~ gsub("_", " ", .x))
species_data_scaled <- sqrt(sqrt(species_data))
species_data_scaled[is.na(species_data_scaled)] <- 0

## Environmental variables
env_data <- biomass_data %>%
  select(Age, Sediment.Rate, PC1, PC2, Pann, TJul)

## Fit RDA model
rda_model <- rda(species_data_scaled ~ Age + Sediment.Rate + PC1 + PC2 + Pann + TJul, data = env_data)
summary(rda_model)

# =====================
# 4. Principal Component Analysis (PCA)
# =====================

## Load plant data
plant_data <- read.csv(paste0(data_dir, "Viridiplantae_PCA_combin_all.csv"))

## Assign lake group
plant_data <- plant_data %>%
  mutate(Lake_group = case_when(
    grepl("Ilirney", Lake, ignore.case = TRUE) ~ "Ilirney",
    grepl("Lele", Lake, ignore.case = TRUE) ~ "Lele",
    grepl("Lama", Lake, ignore.case = TRUE) ~ "Lama",
    grepl("Btoko", Lake, ignore.case = TRUE) ~ "BToko",
    grepl("Ulu", Lake, ignore.case = TRUE) ~ "Ulu",
    grepl("Salmon", Lake, ignore.case = TRUE) ~ "Salmon",
    TRUE ~ NA_character_
  ),
  Lake_group = factor(Lake_group, levels = c("Ilirney", "Lele", "Lama", "BToko", "Salmon", "Ulu")))

## Hellinger transformation and normalization
plant_species <- plant_data[, 3:30]
plant_hellinger <- decostand(plant_species, method = "hellinger")

plant_normalized <- plant_hellinger %>%
  as.data.frame() %>%
  bind_cols(Lake_group = plant_data$Lake_group) %>%
  group_by(Lake_group) %>%
  mutate(across(where(is.numeric), ~ (.-min(.)) / (max(.) - min(.)))) %>%
  ungroup() %>%
  select(-Lake_group)

## PCA using rda (as PCA)
pca_result <- rda(plant_normalized)
summary(pca_result)