# =============================================================================
# Fig. 6 – Comparison of DNA weight and DNA-projected biomass composition
# =============================================================================
#
# Description:
#   This script compares DNA weight percentage with DNA-projected biomass
#   weight percentage across major taxonomic groups and summarizes the
#   relative biomass composition of aquatic and terrestrial taxa.
#
# Panels:
#   Fig. 6a – DNA weight percentage vs. DNA-projected biomass weight percentage
#   Fig. 6b – Relative composition of DNA-projected biomass
#
# Analyses:
#   - Calculation of DNA weight percentage by taxonomic group
#   - Calculation of DNA-projected biomass weight percentage
#   - Comparison of DNA weight and biomass weight distributions
#   - Summary of biomass composition across taxonomic groups
#   - Median and interquartile range visualization
#   - Pseudo-log transformation for low-abundance values
#
# Data:
#   data/Processed_DNA_with_MC_biomass_final.csv
#
# Output:
#   figures
#
# =============================================================================

# ==============================
# Libraries
# ==============================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggthemes)
library(patchwork)
library(readr)

# ==============================
# Paths
# ==============================
data_file <- "data/Processed_DNA_with_MC_biomass_final.csv"
output_dir <- "figures"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================
# Read data
# ==============================
Burial_rate_biomass <- read.csv(data_file)

Burial_rate_biomass[is.na(Burial_rate_biomass)] <- 0


# ==============================
# Add TimeSlice
# ==============================
timeslices <- seq(0, 60, 2)

assign_timeslice <- function(Age) {
  for (i in 1:(length(timeslices) - 1)) {
    if (timeslices[i] <= Age & Age < timeslices[i + 1]) {
      return(paste0(timeslices[i], "-", timeslices[i + 1]))
    }
  }
  
  return(paste0(timeslices[length(timeslices)], "+"))
}

Burial_rate_biomass$timeslice <- sapply(
  Burial_rate_biomass$Age,
  assign_timeslice
)

Burial_rate_biomass$Startage <- as.integer(
  sub("-.*", "", Burial_rate_biomass$timeslice)
)


# ==============================
# Standardize core names
# ==============================
Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    core_name = case_when(
      grepl("Ilirney", Core, ignore.case = TRUE) ~ "Ilirney",
      grepl("Lele", Core, ignore.case = TRUE) ~ "Lele",
      grepl("Lama", Core, ignore.case = TRUE) ~ "Lama",
      grepl("Btoko", Core, ignore.case = TRUE) ~ "BToko",
      grepl("Salmon", Core, ignore.case = TRUE) ~ "Salmon",
      grepl("Ulu", Core, ignore.case = TRUE) ~ "Ulu"
    )
  )


# ==============================
# Taxon order and colors
# ==============================
taxon_order <- c(
  "Aquatic algae",
  "Archaea aquatic",
  "Bacteria aquatic",
  "Fungi aquatic",
  "Metazoa aquatic",
  "Viridiplantae aquatic",
  "Archaea terrestrial",
  "Bacteria terrestrial",
  "Fungi terrestrial",
  "Metazoa terrestrial",
  "Viridiplantae woody",
  "Viridiplantae non woody",
  "Viruses"
)

color_values <- c(
  "#ADD8E6",
  "#00BFFF",
  "#4169E1",
  "#0000FF",
  "#00008B",
  "royalblue4",
  "#98FB98",
  "#50EE50",
  "#32CD32",
  "#228B22",
  "#336433",
  "#006400",
  "#003300"
)


# ============================================================
# Figure 6a: DNA weight percentage vs Biomass weight percentage
# ============================================================

# ------------------------------
# DNA data
# ------------------------------
df_dna <- Burial_rate_biomass %>%
  select(
    "Aquatic_algae_percentage_Start_wt",
    "Bacteria_terrestrial_percentage_Start_wt",
    "Bacteria_aquatic_percentage_Start_wt",
    "Viridiplantae_woody_percentage_Start_wt",
    "Viridiplantae_non_woody_percentage_Start_wt",
    "Viridiplantae_aquatic_percentage_Start_wt",
    "Archaea_terrestrial_percentage_Start_wt",
    "Archaea_aquatic_percentage_Start_wt",
    "Fungi_terrestrial_percentage_Start_wt",
    "Fungi_aquatic_percentage_Start_wt",
    "Metazoa_terrestrial_percentage_Start_wt",
    "Metazoa_aquatic_percentage_Start_wt",
    "Viruses_percentage_Start_wt"
  ) %>%
  mutate(
    across(everything(), ~ . * 100)
  ) %>%
  rename_all(~ gsub("_percentage_Start_wt$", "", .)) %>%
  rename_all(~ gsub("_", " ", .)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Taxon",
    values_to = "Value"
  ) %>%
  filter(Value > 0) %>%
  mutate(
    Taxon = factor(Taxon, levels = taxon_order),
    Source = "DNA weight percentage"
  )


# ------------------------------
# Biomass data
# ------------------------------
df_biomass <- Burial_rate_biomass %>%
  select(
    median_biomass_Aquatic_algae,
    median_biomass_Archaea_aquatic,
    median_biomass_Bacteria_aquatic,
    median_biomass_Fungi_aquatic,
    median_biomass_Metazoa_aquatic,
    median_biomass_Viridiplantae_aquatic,
    median_biomass_Archaea_terrestrial,
    median_biomass_Bacteria_terrestrial,
    median_biomass_Fungi_terrestrial,
    median_biomass_Metazoa_terrestrial,
    median_biomass_Viridiplantae_woody,
    median_biomass_Viridiplantae_non_woody,
    median_biomass_Viruses
  ) %>%
  rename_with(~ gsub("^median_biomass_", "", .)) %>%
  rename_with(~ gsub("_", " ", .)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Taxon",
    values_to = "Value"
  ) %>%
  mutate(
    Value = Value * 100,
    Taxon = factor(Taxon, levels = taxon_order),
    Source = "Biomass weight percentage"
  )


# ------------------------------
# Combine DNA and Biomass
# ------------------------------
df_combined <- bind_rows(
  df_dna,
  df_biomass
)


# ------------------------------
# Summary statistics
# ------------------------------
df_combined_summary <- df_combined %>%
  group_by(Taxon, Source) %>%
  summarise(
    median = median(Value),
    q1 = quantile(Value, 0.25),
    q3 = quantile(Value, 0.75),
    .groups = "drop"
  ) %>%
  mutate(
    Taxon = factor(Taxon, levels = taxon_order)
  )


# ==============================
# Pseudo-log transformation
# ==============================
pseudo_log_trans <- function(base = 10, sigma = 1) {
  
  trans <- function(x) {
    sign(x) * log10(1 + abs(x) / sigma)
  }
  
  inv <- function(x) {
    sign(x) * sigma * (10^abs(x) - 1)
  }
  
  trans_new(
    "pseudo_log",
    trans,
    inv,
    domain = c(0, Inf)
  )
}


# Pseudo-log breaks
breaks <- c(
  1e-6,
  1e-5,
  1e-4,
  1e-3,
  1e-2,
  1e-1
)

color_values_alpha <- alpha(
  color_values,
  0.6
)


# ------------------------------
# DNA vs Biomass weight percentage
# ------------------------------
wt_group <- ggplot(
  df_combined_summary,
  aes(x = Taxon, y = median)
) +
  geom_col(
    aes(
      fill = ifelse(
        Source == "Biomass weight percentage",
        as.character(Taxon),
        NA
      ),
      color = as.character(Taxon),
      group = Source,
      size = ifelse(
        Source == "DNA weight percentage",
        1.2,
        0.5
      )
    ),
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  geom_errorbar(
    aes(
      ymin = q1,
      ymax = q3,
      color = Taxon,
      group = Source
    ),
    position = position_dodge(width = 0.7),
    width = 0.2,
    size = 1.2
  ) +
  scale_fill_manual(
    values = setNames(
      color_values_alpha,
      taxon_order
    ),
    na.value = "transparent"
  ) +
  scale_color_manual(
    values = setNames(
      color_values,
      taxon_order
    ),
    guide = "none"
  ) +
  scale_size_identity() +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-6),
    breaks = breaks,
    labels = c(
      expression(10^{-6}),
      expression(10^{-5}),
      expression(10^{-4}),
      expression(10^{-3}),
      expression(10^{-2}),
      expression(10^{-1})
    )
  ) +
  theme_bw(base_size = 25) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "gray40",
      size = 0.8
    ),
    legend.position = "none"
  ) +
  labs(
    y = "Weight percentage (wt%, log 10 scales)",
    fill = "Taxon",
    color = "Taxon"
  )


# ============================================================
# Figure 6b: Biomass percentage
# ============================================================

df_selected <- Burial_rate_biomass %>%
  select(
    median_biomass_Aquatic_algae_percent,
    median_biomass_Archaea_aquatic_percent,
    median_biomass_Bacteria_aquatic_percent,
    median_biomass_Fungi_aquatic_percent,
    median_biomass_Metazoa_aquatic_percent,
    median_biomass_Viridiplantae_aquatic_percent,
    median_biomass_Archaea_terrestrial_percent,
    median_biomass_Bacteria_terrestrial_percent,
    median_biomass_Fungi_terrestrial_percent,
    median_biomass_Metazoa_terrestrial_percent,
    median_biomass_Viridiplantae_woody_percent,
    median_biomass_Viridiplantae_non_woody_percent,
    median_biomass_Viruses_percent
  ) %>%
  rename_with(~ gsub("^median_biomass_", "", .)) %>%
  rename_with(~ gsub("_percent$", "", .)) %>%
  rename_with(~ gsub("_", " ", .))


# ------------------------------
# Convert to long format
# ------------------------------
df_long <- df_selected %>%
  pivot_longer(
    cols = everything(),
    names_to = "Taxon",
    values_to = "Percentage"
  ) %>%
  mutate(
    Taxon = factor(
      Taxon,
      levels = taxon_order
    )
  )


# ------------------------------
# Summary statistics
# ------------------------------
df_summary <- df_long %>%
  group_by(Taxon) %>%
  summarise(
    median = median(
      Percentage,
      na.rm = TRUE
    ),
    q1 = quantile(
      Percentage,
      0.25,
      na.rm = TRUE
    ),
    q3 = quantile(
      Percentage,
      0.75,
      na.rm = TRUE
    )
  ) %>%
  mutate(
    Taxon = factor(
      Taxon,
      levels = taxon_order
    )
  )


# ------------------------------
# Biomass percentage plot
# ------------------------------
Percentage_wt_group <- ggplot(
  df_summary,
  aes(
    x = Taxon,
    y = median,
    fill = Taxon,
    color = Taxon
  )
) +
  geom_col(
    width = 0.6,
    size = 1.2
  ) +
  geom_errorbar(
    aes(
      ymin = q1,
      ymax = q3
    ),
    width = 0.2,
    size = 1.2
  ) +
  scale_x_discrete(
    limits = taxon_order
  ) +
  coord_cartesian(
    ylim = c(0, 42)
  ) +
  scale_y_continuous(
    name = "Percentage of OC DNA-projected weight percentage (%)"
  ) +
  scale_fill_manual(
    values = alpha(
      color_values,
      0.6
    )
  ) +
  scale_color_manual(
    values = color_values
  ) +
  theme_bw(base_size = 25) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 25,
      angle = 0,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 25
    ),
    axis.title = element_text(
      size = 25
    ),
    axis.title.x = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "gray40",
      size = 0.8
    ),
    axis.ticks.length = unit(
      2,
      "mm"
    ),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )


# ==============================
# Combine plots
# ==============================
combined_plot <- wt_group / Percentage_wt_group

combined_plot <- combined_plot +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(
      size = 25,
      face = "bold"
    )
  )


# ==============================
# Save Figure 6
# ==============================
ggsave(
  file.path(
    output_dir,
    "Fig6_DNA_vs_Biomass.png"
  ),
  combined_plot,
  width = 15,
  height = 20
)

ggsave(
  file.path(
    output_dir,
    "Fig6_DNA_vs_Biomass.pdf"
  ),
  combined_plot,
  width = 15,
  height = 20
)

