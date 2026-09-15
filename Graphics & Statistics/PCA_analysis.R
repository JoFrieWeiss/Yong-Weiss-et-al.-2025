
# ============================================================
# PCA analysis and PCA scores
# ============================================================

library(dplyr)
library(vegan)


# ============================================================
# 1. Load data
# ============================================================

data <- read.csv(
  "data/Comibne_plant_percentage.csv",
  header = TRUE
)


# ============================================================
# 2. Prepare species data
# ============================================================

species_data <- data[, 3:25]

data <- data %>%
  mutate(
    Lake = case_when(
      grepl("Ilirney", Lake, ignore.case = TRUE) ~ "Ilirney",
      grepl("Lele", Lake, ignore.case = TRUE) ~ "Lele",
      grepl("Lama", Lake, ignore.case = TRUE) ~ "Lama",
      grepl("Btoko", Lake, ignore.case = TRUE) ~ "BToko",
      grepl("Ulu", Lake, ignore.case = TRUE) ~ "Ulu",
      grepl("Salmon", Lake, ignore.case = TRUE) ~ "Salmon",
      TRUE ~ NA_character_
    )
  )

lake_levels <- c(
  "Ilirney",
  "Lele",
  "Lama",
  "BToko",
  "Salmon",
  "Ulu"
)

data$Lake_group <- factor(
  data$Lake,
  levels = lake_levels
)


# ============================================================
# 3. Hellinger transformation
# ============================================================

species_data <- decostand(
  species_data,
  method = "hellinger"
)

species_data[is.na(species_data)] <- 0


# ============================================================
# 4. Normalize species data within Lake
# ============================================================

species_data_normalized <- species_data %>%
  as.data.frame() %>%
  bind_cols(
    Lake_group = data$Lake_group
  ) %>%
  group_by(Lake_group) %>%
  mutate(
    across(
      where(is.numeric),
      ~ (.-min(.)) / (max(.) - min(.))
    )
  ) %>%
  ungroup()

species_data_normalized[is.na(species_data_normalized)] <- 0

species_data_normalized <- species_data_normalized %>%
  select(where(is.numeric))


# ============================================================
# 5. PCA
# ============================================================

pca_result <- rda(
  species_data_normalized
)


# ============================================================
# 6. Extract PCA site scores
# ============================================================

pca_scores <- scores(
  pca_result,
  display = "sites"
)

pca_scores <- as.data.frame(pca_scores)


# Add Lake and Age information

pca_scores$Lake <- data$Lake
pca_scores$Age <- data$Age


# ============================================================
# 7. Calculate explained variance
# ============================================================

eigenvalues <- eigenvals(pca_result)

explained_variance <-
  eigenvalues / sum(eigenvalues) * 100

PC1_variance <- explained_variance[1]
PC2_variance <- explained_variance[2]


# ============================================================
# 8. Export PCA scores
# ============================================================

write.csv(
  pca_scores,
  "results/PCA_scores.csv",
  row.names = FALSE
)


# ============================================================
# End of analysis
# ============================================================

