# =============================================================================
# Fig. 2 – OC DNA-projected burial rate
# =============================================================================
#
# Description:
#   This script calculates and visualizes OC DNA-projected burial rate
#   based on DNA-projected biomass.
#
# Panels:
#   Fig. 2  – OC DNA-projected burial rate by lake/core
#   Fig. 2c – OC DNA-projected burial rate across Marine Isotope Stages (MIS)
#
# Data:
#   data/Processed_DNA_with_MC_biomass_final.csv
#
# Output:
#   figures
#
# =============================================================================


# Load required packages

library(tidyverse)

# ============================================================

# 1. Load data

# ============================================================

# Relative path for GitHub reproducibility

Burial_rate_biomass <- read.csv(
  "data/Processed_DNA_with_MC_biomass_final.csv"
)

# Output directory

output_dir <- "figures"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================

# Fig. 2 – Burial rate by lake/core

# ============================================================

# Define Timeslice range

timeslices <- seq(0, 60, 2)

# Assign Age to Timeslice

assign_timeslice <- function(Age) {
  
  for (i in 1:(length(timeslices) - 1)) {
    
    
    if (timeslices[i] <= Age &
        Age < timeslices[i + 1]) {
      
      return(
        paste0(
          timeslices[i],
          "-",
          timeslices[i + 1]
        )
      )
    }
    
    
  }
  
  return(
    paste0(
      timeslices[length(timeslices)],
      "+"
    )
  )
}

# Add Timeslice and Startage

Burial_rate_biomass$timeslice <- sapply(
  Burial_rate_biomass$Age,
  assign_timeslice
)

Burial_rate_biomass$Startage <- as.integer(
  sub(
    "-.*",
    "",
    Burial_rate_biomass$timeslice
  )
)

# ------------------------------------------------------------

# Standardize core names

# ------------------------------------------------------------

Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    core_name = case_when(
      grepl(
        "Ilirney",
        Core.x,
        ignore.case = TRUE
      ) ~ "Ilirney",
      
      
      grepl(
        "Lele",
        Core.x,
        ignore.case = TRUE
      ) ~ "Lele",
      
      grepl(
        "Lama",
        Core.x,
        ignore.case = TRUE
      ) ~ "Lama",
      
      grepl(
        "Btoko",
        Core.x,
        ignore.case = TRUE
      ) ~ "BToko",
      
      grepl(
        "Salmon",
        Core.x,
        ignore.case = TRUE
      ) ~ "Salmon",
      
      grepl(
        "Ulu",
        Core.x,
        ignore.case = TRUE
      ) ~ "Ulu"
    )
    
    
  )

# ------------------------------------------------------------

# Core colours

# ------------------------------------------------------------

core_colors <- c(
  "Ilirney" = "brown",
  "Lele" = "forestgreen",
  "Lama" = "royalblue4",
  "BToko" = "#87CEFF",
  "Salmon" = "darkolivegreen",
  "Ulu" = "orange"
)

# ------------------------------------------------------------

# Prepare data

# ------------------------------------------------------------

Burial_rate_biomass_long <- Burial_rate_biomass %>%
  select(
    Startage,
    BR_Total_median_biomass,
    core_name
  ) %>%
  pivot_longer(
    cols = BR_Total_median_biomass,
    names_to = "Type",
    values_to = "Value"
  )

Burial_rate_biomass_long$core_name <- factor(
  Burial_rate_biomass_long$core_name,
  levels = c(
    "BToko",
    "Ulu",
    "Lama",
    "Ilirney",
    "Salmon",
    "Lele"
  )
)

# ------------------------------------------------------------

# Calculate median and quartiles

# ------------------------------------------------------------

BR_biomass_summary <- Burial_rate_biomass_long %>%
  group_by(core_name) %>%
  summarise(
    median = median(
      Value,
      na.rm = TRUE
    ),
    
    
    Q1 = quantile(
      Value,
      0.25,
      na.rm = TRUE
    ),
    
    Q3 = quantile(
      Value,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
    
  ) %>%
  mutate(
    ymin = Q1,
    ymax = Q3
  )

# ------------------------------------------------------------

# Plot Fig. 2 – Burial rate by lake/core

# ------------------------------------------------------------

bar_BR_biomass_lake <- ggplot(
  BR_biomass_summary,
  aes(
    x = core_name,
    y = median,
    color = core_name,
    fill = core_name
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.7,
    size = 2,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = ymin,
      ymax = ymax
    ),
    width = 0.2,
    size = 2
  ) +
  
  scale_x_discrete(
    name = ""
  ) +
  
  scale_y_continuous(
    name = "OC DNA-projected burial rate\n(g/cm²/year)",
    limits = c(0, 0.008)
  ) +
  
  scale_fill_manual(
    values = alpha(
      core_colors,
      0.6
    )
  ) +
  
  scale_color_manual(
    values = core_colors
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none",
    
    
    axis.text.x = element_text(
      size = 8
    ),
    
    axis.text.y = element_text(
      size = 8
    ),
    
    axis.title = element_text(
      size = 8
    ),
    
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      size = 2
    ),
    
    axis.ticks.length = unit(
      4,
      "mm"
    ),
    
    panel.grid = element_blank(),
    
    panel.background = element_blank()
   
    
  )

# ------------------------------------------------------------

# Save Fig. 2

# ------------------------------------------------------------

ggsave(
  file.path(
    output_dir,
    "Fig2_map.png"
  ),
  bar_BR_biomass_lake,
  width = 5,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "Fig2_map.pdf"
  ),
  bar_BR_biomass_lake,
  width = 5,
  height = 4
)

# ============================================================

# Fig. 2c – MIS 3 / MIS 2 / MIS 1

# ============================================================

Burial_rate_biomass_long <- Burial_rate_biomass %>%
  select(
    Age,
    BR_Total_median_biomass
  ) %>%
  pivot_longer(
    cols = BR_Total_median_biomass,
    names_to = "Type",
    values_to = "Value"
  )

# Assign Marine Isotope Stage

Burial_rate_biomass_long <- Burial_rate_biomass_long %>%
  mutate(
    MIS = case_when(
      Age > 26 ~ "MIS 3 (>26 kyr)",
      Age > 12 & Age <= 26 ~ "MIS 2 (26–12 kyr)",
      Age <= 12 ~ "MIS 1 (12–0 kyr)"
    )
  )

# Set MIS order

Burial_rate_biomass_long$MIS <- factor(
  Burial_rate_biomass_long$MIS,
  levels = c(
    "MIS 3 (>26 kyr)",
    "MIS 2 (26–12 kyr)",
    "MIS 1 (12–0 kyr)"
  )
)

# Calculate median and quartiles

BR_summary <- Burial_rate_biomass_long %>%
  group_by(MIS) %>%
  summarise(
    median = median(
      Value,
      na.rm = TRUE
    ),
    
    
    Q1 = quantile(
      Value,
      0.25,
      na.rm = TRUE
    ),
    
    Q3 = quantile(
      Value,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
    
  ) %>%
  mutate(
    ymin = Q1,
    ymax = Q3
  )

# MIS colours

mis_colors <- c(
  "MIS 3 (>26 kyr)" = "#A6CEE3",
  "MIS 2 (26–12 kyr)" = "#1F78B4",
  "MIS 1 (12–0 kyr)" = "#08306B"
)

# ------------------------------------------------------------

# Plot Fig. 2c

# ------------------------------------------------------------

p_bar <- ggplot(
  BR_summary,
  aes(
    x = MIS,
    y = median,
    fill = MIS,
    color = MIS
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.7,
    size = 1.5,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = ymin,
      ymax = ymax
    ),
    width = 0.2,
    size = 1.2
  ) +
  
  scale_y_continuous(
    name = "OC DNA-projected weight percentage wt%"
  ) +
  
  scale_fill_manual(
    values = mis_colors
  ) +
  
  scale_color_manual(
    values = mis_colors
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none",
    
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title = element_text(
      size = 10
    ),
    
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      size = 1.5
    ),
    
    panel.grid = element_blank()
    
    
  )

# ------------------------------------------------------------

# Save Fig. 2c

# ------------------------------------------------------------

ggsave(
  file.path(
    output_dir,
    "Fig2_map_MIS_barplot.png"
  ),
  p_bar,
  width = 5,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "Fig2_map_MIS_barplot.pdf"
  ),
  p_bar,
  width = 5,
  height = 4
)

# ============================================================

# End of script

# ============================================================
