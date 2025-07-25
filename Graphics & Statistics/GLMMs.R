################################################################################
################################################################################
############################# GLMMs ####################################
################################################################################


# Load required packages
library(lmerTest)
library(dplyr)

# Read the dataset
df <- read.csv("D:/document/DOC/AWI/Lake organic/data/Calculate_all_information_result_unclassfied.csv")

# Select relevant environmental variables
env_df <- df %>%
  select(Number, Age, Lake, Sediment.Rate, PC1, PC2, Pann, TJul)

# Standardize variables within each lake
env_std <- env_df %>%
  group_by(Lake) %>%
  mutate(across(c(Age, Sediment.Rate, PC1, PC2, Pann, TJul), 
                ~ (.-mean(.)) / sd(.), 
                .names = "scaled_{.col}"))

# Extract and log-transform the response variable
dna_df <- df %>%
  select(Number, Age, BR_biomass) %>%
  mutate(BR_Biomass = log10(BR_biomass))

# Merge environmental and response data
combine_data <- merge(dna_df, env_std, by = c("Number", "Age"))

# Fit linear mixed-effects model (Lake as random effect)
fit <- lmer(
  BR_Biomass ~ scaled_Age + scaled_Sediment.Rate + scaled_PC1 + 
    scaled_PC2 + scaled_Pann + scaled_TJul + (1 | Lake), 
  data = combine_data
)

# Output model summary
summary(fit)
