
# =============================================================================
# Fig_4_plot
# Integrated comparison of TOC and DNA-projected biomass
# =============================================================================

# -----------------------------------------------------------------------------
# 0. LIBRARIES
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)


# -----------------------------------------------------------------------------
# 1. PATHS
# -----------------------------------------------------------------------------

PATH_DATA_INPUT <- "data/Processed_DNA_with_MC_biomass.csv"
PATH_OUTPUT_DIR <- "figures/Figure_R1"

if (!dir.exists(PATH_OUTPUT_DIR)) {
  dir.create(PATH_OUTPUT_DIR, recursive = TRUE)
}


# -----------------------------------------------------------------------------
# 2. LOAD DATA
# -----------------------------------------------------------------------------

Burial_rate_biomass <- read.csv(PATH_DATA_INPUT)

Burial_rate_biomass[is.na(Burial_rate_biomass)] <- 0


# -----------------------------------------------------------------------------
# 3. CALCULATE DNA WEIGHT AND BIOMASS
# -----------------------------------------------------------------------------

Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    Total_DNA_weight = rowSums(
      select(., ends_with("_percentage_Start_wt")),
      na.rm = TRUE
    ) * 100
  )


biomass_cols <- c(
  "Aquatic_algae_median_biomass",
  "Archaea_aquatic_median_biomass",
  "Bacteria_aquatic_median_biomass",
  "Fungi_aquatic_median_biomass",
  "Metazoa_aquatic_median_biomass",
  "Viridiplantae_aquatic_median_biomass",
  "Viruses_median_biomass",
  "Archaea_terrestrial_median_biomass",
  "Bacteria_terrestrial_median_biomass",
  "Fungi_terrestrial_median_biomass",
  "Metazoa_terrestrial_median_biomass",
  "Viridiplantae_woody_median_biomass",
  "Viridiplantae_non_woody_median_biomass"
)

aquatic_cols <- biomass_cols[1:6]
terrestrial_cols <- biomass_cols[8:13]


Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    Biomass_Aquatic = rowSums(
      select(., all_of(aquatic_cols)),
      na.rm = TRUE
    ),
    Biomass_Terrestrial = rowSums(
      select(., all_of(terrestrial_cols)),
      na.rm = TRUE
    ),
    Total_median_biomass = rowSums(
      select(., all_of(biomass_cols)),
      na.rm = TRUE
    ),
    Biomass_Aquatic_percentage = ifelse(
      Total_median_biomass > 0,
      Biomass_Aquatic / Total_median_biomass * 100,
      0
    ),
    Biomass_Terrestrial_percentage = ifelse(
      Total_median_biomass > 0,
      Biomass_Terrestrial / Total_median_biomass * 100,
      0
    )
  ) %>%
  mutate(
    across(
      all_of(biomass_cols),
      ~ ifelse(
        Total_median_biomass > 0,
        .x / Total_median_biomass * 100,
        0
      ),
      .names = "{.col}_percent"
    )
  )


# -----------------------------------------------------------------------------
# 4. TIME SLICE AND CORE NAME
# -----------------------------------------------------------------------------

timeslices <- seq(0, 60, 2)

assign_timeslice <- function(Age) {
  for (i in 1:(length(timeslices) - 1)) {
    if (timeslices[i] <= Age & Age < timeslices[i + 1]) {
      return(paste0(timeslices[i], "-", timeslices[i + 1]))
    }
  }
  
  return(paste0(timeslices[length(timeslices)], "+"))
}

Burial_rate_biomass$Age <- as.numeric(
  as.character(Burial_rate_biomass$Age)
)

Burial_rate_biomass$timeslice <- sapply(
  Burial_rate_biomass$Age,
  assign_timeslice
)

Burial_rate_biomass$Startage <- as.integer(
  sub("-.*", "", Burial_rate_biomass$timeslice)
)


Burial_rate_biomass <- Burial_rate_biomass %>%
  mutate(
    core_name = case_when(
      grepl("Ilirney", Core.x, ignore.case = TRUE) ~ "Ilirney",
      grepl("Lele", Core.x, ignore.case = TRUE) ~ "Lele",
      grepl("Lama", Core.x, ignore.case = TRUE) ~ "Lama",
      grepl("Btoko", Core.x, ignore.case = TRUE) ~ "BToko",
      grepl("Salmon", Core.x, ignore.case = TRUE) ~ "Salmon",
      grepl("Ulu", Core.x, ignore.case = TRUE) ~ "Ulu",
      TRUE ~ NA_character_
    )
  )


# -----------------------------------------------------------------------------
# 5. BASE THEME
# -----------------------------------------------------------------------------

base_theme <- theme_bw(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "gray50",
      size = 0.7
    )
  )


# -----------------------------------------------------------------------------
# 6. CORE-WISE COMPARISON PLOTS
# -----------------------------------------------------------------------------

plot_list_cores <- list()

core_list <- unique(Burial_rate_biomass$Number)

for (core_id in core_list) {
  
  core_data <- Burial_rate_biomass %>%
    filter(Number == core_id)
  
  if (nrow(core_data) == 0) {
    next
  }
  
  lake_name <- core_data$Lake[1]
  
  df_final <- core_data %>%
    mutate(
      C_Factor = case_when(
        grepl("Bacteria", Type) ~ 0.52,
        grepl("Viruses", Type) ~ 0.68,
        grepl("Archaea", Type) ~ 0.50,
        grepl("woody", Type) ~ 0.48,
        TRUE ~ 0.50
      ),
      Age_Group = round(Age, 4)
    ) %>%
    group_by(Age_Group) %>%
    summarise(
      Age = mean(Age),
      Sim_Med = sum(
        median_biomass * C_Factor,
        na.rm = TRUE
      ) * 100,
      Sim_Q1 = sum(
        q1_biomass * C_Factor,
        na.rm = TRUE
      ) * 100,
      Sim_Q3 = sum(
        q3_biomass * C_Factor,
        na.rm = TRUE
      ) * 100
    ) %>%
    ungroup()
  
  
  df_meas <- core_data %>%
    select(Age, TOC) %>%
    rename(Measured_TOC = TOC) %>%
    filter(!is.na(Measured_TOC)) %>%
    mutate(
      Meas_Min = Measured_TOC * 0.95,
      Meas_Max = Measured_TOC * 1.05
    )
  
  
  if (nrow(inner_join(df_final, df_meas, by = "Age")) > 2) {
    
    df_test <- inner_join(
      df_final,
      df_meas,
      by = "Age"
    )
    
    corr_res <- cor.test(
      df_test$Sim_Med,
      df_test$Measured_TOC,
      method = "spearman"
    )
    
    wilcox_res <- wilcox.test(
      df_test$Sim_Med,
      df_test$Measured_TOC,
      paired = TRUE
    )
    
    median_ratio <- median(
      df_test$Sim_Med / df_test$Measured_TOC,
      na.rm = TRUE
    )
    
    mag_label <- ifelse(
      wilcox_res$p.value > 0.05,
      "Consistent",
      "Different"
    )
    
    subtitle_text <- paste0(
      "Trend: R=",
      round(corr_res$estimate, 2),
      " (p=",
      format.pval(corr_res$p.value, digits = 2),
      ") | Mag: Ratio ",
      round(median_ratio, 2),
      "x (p=",
      format.pval(wilcox_res$p.value, digits = 2),
      " ",
      mag_label,
      ")"
    )
    
  } else {
    
    subtitle_text <- "Insufficient data"
    
  }
  
  
  p_corr <- ggplot() +
    geom_ribbon(
      data = df_meas,
      aes(
        x = Age,
        ymin = Meas_Min,
        ymax = Meas_Max,
        fill = "Measured TOC (±5%)"
      ),
      alpha = 0.25
    ) +
    geom_smooth(
      data = df_meas,
      aes(
        x = Age,
        y = Measured_TOC,
        color = "Measured TOC (±5%)"
      ),
      method = "loess",
      span = 0.5,
      se = FALSE
    ) +
    geom_ribbon(
      data = df_final,
      aes(
        x = Age,
        ymin = Sim_Q1,
        ymax = Sim_Q3,
        fill = "Simulated Biomass (IQR)"
      ),
      alpha = 0.15
    ) +
    geom_smooth(
      data = df_final,
      aes(
        x = Age,
        y = Sim_Med,
        color = "Simulated Biomass (IQR)"
      ),
      method = "loess",
      span = 0.5,
      se = FALSE
    ) +
    scale_color_manual(
      name = "",
      values = c(
        "Measured TOC (±5%)" = "#5d2c04",
        "Simulated Biomass (IQR)" = "#BC8F8F"
      )
    ) +
    scale_fill_manual(
      name = "",
      values = c(
        "Measured TOC (±5%)" = "#5d2c04",
        "Simulated Biomass (IQR)" = "#BC8F8F"
      )
    ) +
    annotation_logticks(sides = "l") +
    theme_bw(base_size = 8) +
    theme(
      legend.position = "top",
      panel.grid = element_blank()
    ) +
    labs(
      title = paste("", lake_name),
      subtitle = subtitle_text,
      y = "Weight percentage (wt%)",
      x = "Age (k yrs)"
    )
  
  plot_list_cores[[as.character(core_id)]] <- p_corr
}


# -----------------------------------------------------------------------------
# 7. CORE GRID
# -----------------------------------------------------------------------------

final_grid <- wrap_plots(
  plot_list_cores,
  ncol = 3
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Comparison of Simulated Biomass vs. Measured TOC",
    subtitle = "Trend: Pearson Correlation | Magnitude: Median Ratio & Wilcoxon Test",
    theme = theme(
      legend.position = "top",
      legend.direction = "horizontal",
      plot.title = element_text(
        size = 18,
        face = "bold"
      )
    )
  )


# -----------------------------------------------------------------------------
# 8. TOC VS DNA-PROJECTED BIOMASS BAR PLOT
# -----------------------------------------------------------------------------

Burial_rate_biomass_long_total <- Burial_rate_biomass %>%
  select(TOC, Total_median_biomass) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(
    Type = ifelse(
      Type == "Total_median_biomass",
      "OC DNA-projected",
      Type
    ),
    Type = factor(
      Type,
      levels = c("TOC", "OC DNA-projected")
    ),
    Value = ifelse(
      Type == "OC DNA-projected",
      Value * 100,
      Value
    )
  )


TOC_biomass_total_bar <- Burial_rate_biomass_long_total %>%
  group_by(Type) %>%
  summarise(
    median_val = median(Value),
    q1 = quantile(Value, 0.25),
    q3 = quantile(Value, 0.75),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = Type,
      y = median_val,
      fill = Type,
      color = Type
    )
  ) +
  geom_bar(
    stat = "identity",
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    width = 0.2,
    size = 1.2,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "TOC" = "#5d2c04",
      "OC DNA-projected" = "#BC8F8F"
    )
  ) +
  scale_color_manual(
    values = c(
      "TOC" = "#5d2c04",
      "OC DNA-projected" = "#BC8F8F"
    )
  ) +
  labs(
    y = "Weight percentage (%)",
    x = NULL,
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 4)) +
  base_theme +
  theme(legend.position = "none")


# -----------------------------------------------------------------------------
# 9. FINAL FIGURE 4.1
# -----------------------------------------------------------------------------

final_4panel <- (
  TOC_biomass_total_bar / final_grid
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Integrated Comparison of TOC and DNA-Projected Biomass",
    subtitle = "Core-wise trends and habitat partitioning",
    tag_levels = "a",
    theme = theme(
      legend.position = "bottom",
      plot.title = element_text(
        size = 16,
        face = "bold"
      ),
      plot.subtitle = element_text(size = 11)
    )
  )


# -----------------------------------------------------------------------------
# 10. SAVE
# -----------------------------------------------------------------------------

ggsave(
  file.path(PATH_OUTPUT_DIR, "Fig_4_plot.png"),
  final_4panel,
  width = 36,
  height = 30,
  units = "cm",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(PATH_OUTPUT_DIR, "Fig_4_plot.pdf"),
  final_4panel,
  width = 36,
  height = 30,
  units = "cm",
  dpi = 300,
  bg = "white"
)

