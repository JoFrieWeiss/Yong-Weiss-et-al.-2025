# ============================================================
# Burial Rate and Biomass Analysis
# ============================================================
#
# This script performs the statistical analyses of:
#
# 1. Total biomass burial rate
# 2. Aquatic biomass percentage
#
# The processed dataset already contains all calculated
# biomass, percentage, environmental, and burial-rate variables.
#
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(dplyr)
library(glmmTMB)
library(ggeffects)


# ============================================================
# 2. Load final processed dataset
# ============================================================

Burial_rate_biomass <- read.csv(
  "data/Processed_DNA_with_MC_biomass_final.csv",
  stringsAsFactors = FALSE
)


# ============================================================
# 3. Check the dataset
# ============================================================

str(Burial_rate_biomass)

dim(Burial_rate_biomass)

summary(Burial_rate_biomass)

colnames(Burial_rate_biomass)

# Number of observations per Lake
table(Burial_rate_biomass$Lake)

# Check missing values
colSums(is.na(Burial_rate_biomass))


# ============================================================
# 4. Standardize environmental predictors within Lake
# ============================================================
#
# Continuous environmental predictors are standardized
# separately within each Lake.
#
# Standardization:
# (value - mean) / standard deviation
#
# ============================================================

combine_data <- Burial_rate_biomass %>%
  group_by(Lake) %>%
  mutate(
    
    scaled_Age =
      as.numeric(scale(Age)),
    
    scaled_Sediment.Rate =
      as.numeric(scale(Sediment.Rate)),
    
    scaled_PC1 =
      as.numeric(scale(PC1)),
    
    scaled_PC2 =
      as.numeric(scale(PC2)),
    
    scaled_Pann =
      as.numeric(scale(Pann)),
    
    scaled_TJul =
      as.numeric(scale(TJul)),
    
    scaled_Average.read.length.of.sample =
      as.numeric(
        scale(Average.read.length.of.sample)
      )
    
  ) %>%
  ungroup()


# ============================================================
# 5. TOTAL BIOMASS BURIAL-RATE MODEL
# ============================================================

fit_BR <- glmmTMB(
  
  BR_Total_median_biomass ~
    scaled_Age +
    scaled_Sediment.Rate +
    scaled_PC1 +
    scaled_PC2 +
    scaled_Pann +
    scaled_TJul +
    scaled_Average.read.length.of.sample +
    (1 | Lake),
  
  family = Gamma(link = "log"),
  
  data = combine_data
)


# ============================================================
# 6. Total biomass model summary
# ============================================================

summary(fit_BR)


# ============================================================
# 7. AQUATIC BIOMASS PERCENTAGE MODEL
# ============================================================

# Check the response variable

summary(
  combine_data$Biomass_Aquatic_percentage
)

range(
  combine_data$Biomass_Aquatic_percentage,
  na.rm = TRUE
)


# ------------------------------------------------------------
# Aquatic biomass percentage GLMM
# ------------------------------------------------------------
#
# This follows the model structure used in the original
# analysis.
#
# No family is specified here because the original model
# did not specify one.
#
# ------------------------------------------------------------

fit_aquatic <- glmmTMB(
  
  Biomass_Aquatic_percentage ~
    scaled_Age +
    scaled_Sediment.Rate +
    scaled_PC1 +
    scaled_PC2 +
    scaled_Pann +
    scaled_TJul +
    scaled_Average.read.length.of.sample +
    (1 | Lake),
  
  data = combine_data
)


# ============================================================
# 8. Aquatic biomass percentage model summary
# ============================================================

summary(fit_aquatic)


# ============================================================
# END OF ANALYSIS
# ============================================================