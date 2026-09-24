# =============================================================================
# Fig. 7 – OC DNA-projected biomass weight percentage
# =============================================================================
#
# Description:
#   This script visualizes the temporal variation in OC DNA-projected
#   biomass weight percentage across major taxonomic groups during
#   the last 30 kyr.
#
# Panels:
#   - Age-dependent OC DNA-projected biomass weight percentage
#   - Taxon-specific temporal patterns with independent y-axis scales
#
# Analyses:
#   - Filtering samples younger than 30 kyr
#   - Assignment of samples to 2-kyr age slices
#   - Standardization of lake core names
#   - Conversion of biomass data to long format
#   - Calculation of median and interquartile range for each age slice
#   - Taxon-specific visualization of OC DNA-projected biomass percentage
#
# Data:
#   data/Processed_DNA_with_MC_biomass_final.csv
#
# Output:
#   figures
#
# =============================================================================

# ---------------------------
# 1. Load packages
# ---------------------------

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(scales)


# ---------------------------
# 2. File paths
# ---------------------------

PATH_DATA_INPUT <- "data/Processed_DNA_with_MC_biomass_final.csv"

PATH_FIGURE_PNG <- "figures/Fig7_percent_wt.png"
PATH_FIGURE_PDF <- "figures/Fig7_percent_wt.pdf"


# ---------------------------
# 3. Load data
# ---------------------------

Burial_rate_biomass <- read.csv(
  PATH_DATA_INPUT,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ---------------------------
# 4. Filter data
# ---------------------------

Burial_rate_biomass <- Burial_rate_biomass %>%
  filter(Age < 30)
# ---------------------------

# 3. Define age timeslices

# ---------------------------

timeslices <- seq(0, 30, by = 2)

# Explicitly include the upper boundary

timeslices <- c(timeslices, 30)

# Assign each sample to a 2-kyr age bin

assign_timeslice <- function(Age) {
  
  for (i in seq_len(length(timeslices) - 1)) {
    
    
    if (timeslices[i] <= Age && Age < timeslices[i + 1]) {
      return(
        paste0(timeslices[i], "-", timeslices[i + 1])
      )
    }
    
    
  }
  
  return(
    paste0(timeslices[length(timeslices)], "+")
  )
}

# Apply age-bin assignment

Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    timeslice = sapply(Age, assign_timeslice),
    Startage = as.integer(sub("-.*", "", timeslice))
  )

# ---------------------------

# 4. Standardize core names

# ---------------------------

Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    core_name = case_when(
      grepl("Ilirney", Core.x, ignore.case = TRUE) ~ "Illirney",
      grepl("Lele",    Core.x, ignore.case = TRUE) ~ "Lele",
      grepl("Lama",    Core.x, ignore.case = TRUE) ~ "Lama",
      grepl("Btoko",   Core.x, ignore.case = TRUE) ~ "BToko",
      grepl("Salmon",  Core.x, ignore.case = TRUE) ~ "Salmon",
      grepl("Ulu",     Core.x, ignore.case = TRUE) ~ "Ulu",
      TRUE ~ NA_character_
    )
  )

# ---------------------------

# 5. Define taxonomic order

# ---------------------------

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

# ---------------------------

# 6. Define colors

# ---------------------------

color_values <- c(
  "#ADD8E6",   # Aquatic algae
  "#00BFFF",   # Archaea aquatic
  "#4169E1",   # Bacteria aquatic
  "#0000FF",   # Fungi aquatic
  "#00008B",   # Metazoa aquatic
  "royalblue4",# Viridiplantae aquatic
  "#98FB98",   # Archaea terrestrial
  "#50EE50",   # Bacteria terrestrial
  "#32CD32",   # Fungi terrestrial
  "#228B22",   # Metazoa terrestrial
  "#336433",   # Viridiplantae woody
  "#006400",   # Viridiplantae non woody
  "#003300"    # Viruses
)

color_values_alpha <- alpha(color_values, 0.6)

# ---------------------------

# 7. Convert data to long format

# ---------------------------

biomass_columns <- c(
  "Aquatic_algae_median_biomass_percent",
  "Archaea_aquatic_median_biomass_percent",
  "Bacteria_aquatic_median_biomass_percent",
  "Fungi_aquatic_median_biomass_percent",
  "Metazoa_aquatic_median_biomass_percent",
  "Viridiplantae_aquatic_median_biomass_percent",
  "Viruses_median_biomass_percent",
  "Archaea_terrestrial_median_biomass_percent",
  "Bacteria_terrestrial_median_biomass_percent",
  "Fungi_terrestrial_median_biomass_percent",
  "Metazoa_terrestrial_median_biomass_percent",
  "Viridiplantae_woody_median_biomass_percent",
  "Viridiplantae_non_woody_median_biomass_percent"
)

Burial_rate_biomass_long <- Burial_rate_biomass %>%
  pivot_longer(
    cols = all_of(biomass_columns),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(
    Type = str_remove(Type, "median_biomass_"),
    Type = str_remove(Type, "*percent$"),
    Type = str_replace_all(Type, "*", " "),
    Type = factor(Type, levels = taxon_order)
  )

# ---------------------------

# 8. Plot function

# ---------------------------

plot_median_line_biomass <- function(
    data,
    y_label = "OC DNA projected weight percentage (wt%)",
    show_legend = FALSE
) {
  
  # Calculate Q1, median and Q3 for each age slice
  
  summary_data <- data %>%
    group_by(Startage, Type) %>%
    summarise(
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Startage = factor(
        Startage,
        levels = sort(unique(Startage))
      )
    )
  
  # Generate figure
  
  ggplot(summary_data) +
    
    
  geom_bar(
    aes(
      x = Startage,
      y = Median,
      fill = Type,
      color = Type
    ),
    stat = "identity",
    width = 0.6,
    position = position_dodge(width = 0.001),
    linewidth = 0.3
  ) +
    
    geom_errorbar(
      aes(
        x = Startage,
        ymin = Q1,
        ymax = Q3,
        color = Type
      ),
      width = 0.2,
      linewidth = 1.2,
      position = position_dodge(width = 0.001)
    ) +
    
    facet_wrap(
      ~Type,
      scales = "free_y",
      ncol = 1
    ) +
    
    scale_x_discrete(
      name = "Age (Kyr)"
    ) +
    
    scale_y_continuous(
      name = y_label,
      sec.axis = sec_axis(
        ~ .,
        name = y_label
      )
    ) +
    
    scale_fill_manual(
      values = color_values_alpha,
      breaks = taxon_order
    ) +
    
    scale_color_manual(
      values = color_values,
      breaks = taxon_order
    ) +
    
    theme_bw() +
    
    theme(
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.title = element_text(size = 25),
      
      strip.text = element_text(size = 10),
      strip.background = element_blank(),
      strip.text.position = "top",
      
      panel.border = element_blank(),
      panel.grid = element_blank(),
      
      axis.ticks.length = unit(2, "mm"),
      
      legend.position = ifelse(
        show_legend,
        "right",
        "none"
      ),
      
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black")
    )
  
  
}

# ---------------------------

# 9. Generate Figure 7

# ---------------------------

fig7 <- plot_median_line_biomass(
  Burial_rate_biomass_long
)

# ---------------------------

# 10. Export figure

# ---------------------------

# PNG

ggsave(
  filename = "Fig7_percent_wt.png",
  plot = fig7,
  width = 15,
  height = 20,
  units = "in",
  dpi = 300
)

# PDF

ggsave(
  filename = "Fig7_percent_wt.pdf",
  plot = fig7,
  width = 15,
  height = 20,
  units = "in"
)
