# =============================================================================
# Fig. 3 – DNA weight, biomass, TOC, burial rate, biomass sources,
#           July temperature, and PC1
# =============================================================================
#
# Description:
#   This script analyzes and visualizes temporal variations in DNA weight,
#   DNA-projected biomass, TOC, OC DNA-projected burial rate, aquatic and
#   terrestrial biomass, July temperature, and PC1.
#
# Panels:
#   Fig. 3a – DNA weight
#   Fig. 3b – DNA-projected biomass
#   Fig. 3c – TOC
#   Fig. 3d – OC DNA-projected burial rate
#   Fig. 3e – Aquatic and terrestrial biomass
#   Fig. 3f – July temperature (TJul)
#   Fig. 3g – PC1
#
# Data:
#   data/Processed_DNA_with_MC_biomass_final.csv
#
# Output:
#   figures
#
# =============================================================================


library(tidyverse)
library(corit)
library(zoo)
library(cowplot)
library(patchwork)
library(scales)

# ------------------------------------------------------------

# 2. Load data

# ------------------------------------------------------------

# Relative path for GitHub reproducibility

Burial_rate_biomass <- read.csv(
  "data/Processed_DNA_with_MC_biomass_final.csv"
)

# Output directory

output_dir <- "figures"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

# ============================================================

# Fig. 3

# ============================================================

# ------------------------------------------------------------

# 3. Define Timeslices

# ------------------------------------------------------------

timeslices <- seq(0, 30, 2)
timeslices <- c(timeslices, 30)

assign_timeslice <- function(Age) {
  
  for (i in 1:(length(timeslices) - 1)) {
    
    
    if (
      timeslices[i] <= Age &
      Age < timeslices[i + 1]
    ) {
      
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

# 4. Standardize core names

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
      ) ~ "Ulu",
      
      TRUE ~ as.character(Core.x)
    )
    
    
  )

# ============================================================

# 5. Generate main Fig. 3 panels

# ============================================================

# ------------------------------------------------------------

# 5.1 General plotting function

# ------------------------------------------------------------

generate_plot <- function(
    data,
    column_name,
    plot_title,
    y_label,
    y_limit,
    color_code
) {
  
  data_long <- data %>%
    select(
      Startage,
      all_of(column_name),
      core_name
    ) %>%
    pivot_longer(
      cols = all_of(column_name),
      names_to = "Type",
      values_to = "Value"
    )
  
  # Convert to percentage
  
  if (
    column_name %in%
    c(
      "Total_median_biomass",
      "Total_DNA_weight"
    )
  ) {
    
    
    data_long$Value <-
      data_long$Value * 100
    
    
  }
  
  # Calculate median and IQR
  
  summary_data <- data_long %>%
    group_by(
      Startage,
      Type
    ) %>%
    summarise(
      
      
      median = median(
        Value,
        na.rm = TRUE
      ),
      
      q1 = quantile(
        Value,
        0.25,
        na.rm = TRUE
      ),
      
      q3 = quantile(
        Value,
        0.75,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      ymin = q1,
      ymax = q3
    )
  
  
  # Plot
  
  ggplot(
    summary_data,
    aes(
      x = factor(Startage),
      y = median,
      fill = Type
    )
  ) +
    
    
  geom_col(
    aes(
      color = Type
    ),
    position = position_dodge(0.8),
    width = 0.7,
    size = 0.1
  ) +
    
    geom_errorbar(
      aes(
        ymin = ymin,
        ymax = ymax,
        colour = Type
      ),
      position = position_dodge(0.8),
      width = 0.25,
      size = 1.2
    ) +
    
    scale_y_continuous(
      name = y_label,
      position = "right"
    ) +
    
    scale_x_discrete(
      name = "Age (Kyr)"
    ) +
    
    coord_cartesian(
      ylim = c(
        min(
          summary_data$ymin,
          na.rm = TRUE
        ),
        y_limit
      )
    ) +
    
    scale_fill_manual(
      values = alpha(
        color_code,
        0.6
      )
    ) +
    
    scale_color_manual(
      values = color_code
    ) +
    
    theme_bw() +
    
    theme(
      
      panel.grid = element_blank(),
      
      legend.position = "none",
      
      legend.title = element_blank(),
      
      axis.text.x = element_blank(),
      
      axis.title.x = element_blank(),
      
      axis.ticks.x = element_blank(),
      
      axis.text.y = element_text(
        size = 8
      ),
      
      axis.title = element_text(
        size = 8
      ),
      
      panel.border = element_rect(
        fill = NA,
        colour = "grey",
        size = 0.8
      ),
      
      axis.ticks.length = unit(
        2,
        "mm"
      ),
      
      axis.ticks = element_line(
        size = 1
      )
    )
  
  
}

# ------------------------------------------------------------

# 5.2 Plot parameters

# ------------------------------------------------------------

plot_params <- list(
  
  list(
    column_name = "Total_DNA_weight",
    plot_title = "DNA Weight",
    y_label = "DNA weight percentage\n(wt%)",
    y_limit = 0.008,
    color_code = "#DEB887"
  ),
  
  list(
    column_name = "Total_median_biomass",
    plot_title = "Biomass Weight",
    y_label = "OC DNA-projected\n(wt%)",
    y_limit = 10,
    color_code = "#BC8F8F"
  ),
  
  list(
    column_name = "TOC",
    plot_title = "TOC",
    y_label = "TOC (wt%)",
    y_limit = 10,
    color_code = "#5d2c04"
  ),
  
  list(
    column_name = "BR_Total_median_biomass",
    plot_title = "Biomass Burial Rate",
    y_label = "OC DNA-projected burial\nrate (g/cm²/years)",
    y_limit = 0.004,
    color_code = "#2c1503"
  )
)

# ------------------------------------------------------------

# 5.3 Generate plots

# ------------------------------------------------------------

plot_list <- lapply(
  plot_params,
  function(params) {
    
    
    generate_plot(
      Burial_rate_biomass,
      params$column_name,
      params$plot_title,
      params$y_label,
      params$y_limit,
      params$color_code
    )
    
    
  }
)

# Extract individual plots

DNA_plot <- plot_list[[1]]

Biomass_plot <- plot_list[[2]]

TOC_plot <- plot_list[[3]]

BR_Biomass_plot <- plot_list[[4]]

# ============================================================

# 6. Aquatic vs. Terrestrial biomass

# ============================================================

prepare_burial_data <- function(data) {
  
  data %>%
    
    
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
      ) ~ "Ulu",
      
      TRUE ~ NA_character_
    )
  ) %>%
    
    select(
      Age,
      Biomass_Aquatic_percentage,
      Biomass_Terrestrial_percentage,
      core_name,
      Startage,
      timeslice
    ) %>%
    
    pivot_longer(
      cols = c(
        Biomass_Aquatic_percentage,
        Biomass_Terrestrial_percentage
      ),
      names_to = "Type",
      values_to = "Value"
    ) %>%
    
    mutate(
      Type = recode(
        Type,
        "Biomass_Aquatic_percentage" = "Aquatic",
        "Biomass_Terrestrial_percentage" = "Terrestrial"
      ),
      
      Type = factor(
        Type,
        levels = c(
          "Aquatic",
          "Terrestrial"
        )
      )
    )
  
  
}

prepared_data <- prepare_burial_data(
  Burial_rate_biomass
)

# ------------------------------------------------------------

# Aquatic / Terrestrial colours

# ------------------------------------------------------------

type_colors <- c(
  "Aquatic" = "#1f78b4",
  "Terrestrial" = "#33a02c"
)

type_colors_alpha <- alpha(
  type_colors,
  0.6
)

# ------------------------------------------------------------

# Calculate summary statistics

# ------------------------------------------------------------

summary_biomass <- prepared_data %>%
  
  group_by(
    Startage,
    Type
  ) %>%
  
  summarise(
    
    
    median = median(
      Value,
      na.rm = TRUE
    ),
    
    q1 = quantile(
      Value,
      0.25,
      na.rm = TRUE
    ),
    
    q3 = quantile(
      Value,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
    
  ) %>%
  
  mutate(
    ymin = q1,
    ymax = q3
  )

# ------------------------------------------------------------

# Aquatic vs. Terrestrial plot

# ------------------------------------------------------------

ter_aqu_biomass <- ggplot(
  summary_biomass,
  aes(
    x = factor(Startage),
    y = median,
    fill = Type
  )
) +
  
  geom_col(
    aes(
      color = Type
    ),
    position = position_dodge(0.8),
    width = 0.7,
    size = 0.1
  ) +
  
  geom_errorbar(
    aes(
      ymin = ymin,
      ymax = ymax,
      colour = Type
    ),
    position = position_dodge(0.8),
    width = 0.25,
    size = 1.2
  ) +
  
  scale_fill_manual(
    values = type_colors_alpha,
    name = NULL
  ) +
  
  scale_color_manual(
    values = type_colors
  ) +
  
  scale_y_continuous(
    name = "OC DNA-projected\nweight percentage (wt%)",
    position = "right",
    limits = c(0, 100)
  ) +
  
  scale_x_discrete(
    name = "Age (Kyr)"
  ) +
  
  theme_bw() +
  
  theme(
    
    
    legend.position = c(1, 0),
    
    legend.justification = c(
      1,
      0
    ),
    
    legend.background = element_rect(
      fill = "white",
      color = "black"
    ),
    
    legend.text = element_text(
      size = 20
    ),
    
    legend.key.size = unit(
      0.8,
      "cm"
    ),
    
    axis.text.x = element_blank(),
    
    axis.title.x = element_blank(),
    
    axis.ticks.x = element_blank(),
    
    panel.border = element_rect(
      fill = NA,
      colour = "grey",
      size = 0.8
    ),
    
    axis.ticks.length = unit(
      2,
      "mm"
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      size = 1
    )
    
    
  )

# ============================================================

# 7. Interpolation of TJul

# ============================================================

# Lakes included in the analysis

selected_cores <- unique(
  Burial_rate_biomass$Lake
)

# ------------------------------------------------------------

# Interpolate TJul for each lake

# ------------------------------------------------------------

interpola_cores_filtered_corit <- NULL

for (i in 1:length(selected_cores)) {
  
  print(
    paste0(
      i,
      "/",
      length(selected_cores)
    )
  )
  
  selected_cores_sub <-
    Burial_rate_biomass %>%
    filter(
      Lake == selected_cores[i]
    ) %>%
    mutate(
      Age1 = Age * 1000
    )
  
  selected_cores_sub1 <- zoo(
    selected_cores_sub$TJul,
    order.by = selected_cores_sub$Age1
  )
  
  selected_cores_sub_after <- InterpolationMethod(
    selected_cores_sub1,
    fc = 1 / 100,
    dt = 100,
    tail(
      selected_cores_sub$Age1,
      1
    ),
    int.method = "linear",
    appliedFilter = "gauss",
    k = 5
  )
  
  selected_cores_sub_df <- data.frame(
    Age = time(
      selected_cores_sub_after
    ),
    value = as.vector(
      selected_cores_sub_after
    )
  )
  
  timesteps <- seq(
    0,
    30000,
    1000
  )
  
  selected_cores_sub_ts <-
    selected_cores_sub_df %>%
    filter(
      Age %in% timesteps
    ) %>%
    drop_na()
  
  interpola_df <- data.frame(
    lake = unique(
      selected_cores_sub$Lake
    ),
    selected_cores_sub_ts
  )
  
  interpola_cores_filtered_corit <-
    rbind(
      interpola_cores_filtered_corit,
      interpola_df
    )
}

# ------------------------------------------------------------

# Calculate median TJul

# ------------------------------------------------------------

median_value_TJul <-
  interpola_cores_filtered_corit %>%
  
  arrange(Age) %>%
  
  group_by(Age) %>%
  
  mutate(
    median_TJul = median(value)
  ) %>%
  
  ungroup() %>%
  
  select(
    Age,
    median_TJul
  ) %>%
  
  distinct()

# ------------------------------------------------------------

# Regression and Age = 0 prediction

# ------------------------------------------------------------

median_values <- median_value_TJul %>%
  mutate(
    Age = Age / 1000
  )

lm_TJul <- lm(
  median_TJul ~ Age,
  data = median_values
)

predicted_TJul <- predict(
  lm_TJul,
  newdata = data.frame(
    Age = 0
  )
)

predicted_row <- data.frame(
  Age = 0,
  median_TJul = predicted_TJul
)

median_values <- rbind(
  predicted_row,
  median_values
) %>%
  arrange(Age)

# ------------------------------------------------------------

# TJul plot

# ------------------------------------------------------------

TJul_plot <- ggplot(
  median_values,
  aes(
    x = Age,
    y = median_TJul
  )
) +
  
  geom_line(
    color = "#c5272d",
    size = 1.2
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      28,
      2
    ),
    limits = c(
      0,
      28
    )
  ) +
  
  coord_cartesian(
    ylim = c(
      9,
      13
    )
  ) +
  
  scale_y_continuous(
    name = "July temperature (°C)",
    position = "right"
  ) +
  
  theme_minimal() +
  
  theme(
    
    
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    panel.background = element_rect(
      fill = "white",
      color = "white"
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = "white"
    ),
    
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
      colour = "grey",
      size = 0.8
    ),
    
    axis.line = element_line(
      color = "grey",
      size = 0.8
    ),
    
    axis.ticks.length = unit(
      2,
      "mm"
    ),
    
    axis.ticks = element_line(
      color = "black",
      size = 1
    ),
    
    text = element_text(
      size = 25
    )
    
    
  )

# ============================================================

# 8. Interpolation of PC1

# ============================================================

interpola_cores_filtered_corit_PC1 <- NULL

for (i in 1:length(selected_cores)) {
  
  print(
    paste0(
      "PC1 ",
      i,
      "/",
      length(selected_cores)
    )
  )
  
  selected_cores_sub <-
    Burial_rate_biomass %>%
    filter(
      Lake == selected_cores[i]
    ) %>%
    mutate(
      Age1 = Age * 1000
    )
  
  selected_cores_sub1 <- zoo(
    selected_cores_sub$PC1,
    order.by = selected_cores_sub$Age1
  )
  
  selected_cores_sub_after <- InterpolationMethod(
    selected_cores_sub1,
    fc = 1 / 100,
    dt = 100,
    tail(
      selected_cores_sub$Age1,
      1
    ),
    int.method = "linear",
    appliedFilter = "gauss",
    k = 5
  )
  
  selected_cores_sub_df <- data.frame(
    Age = time(
      selected_cores_sub_after
    ),
    value = as.vector(
      selected_cores_sub_after
    )
  )
  
  timesteps <- seq(
    0,
    30000,
    1000
  )
  
  selected_cores_sub_ts <-
    selected_cores_sub_df %>%
    filter(
      Age %in% timesteps
    ) %>%
    drop_na()
  
  interpola_df <- data.frame(
    lake = unique(
      selected_cores_sub$Lake
    ),
    selected_cores_sub_ts
  )
  
  interpola_cores_filtered_corit_PC1 <-
    rbind(
      interpola_cores_filtered_corit_PC1,
      interpola_df
    )
}

# ------------------------------------------------------------

# Calculate median PC1

# ------------------------------------------------------------

median_value_PC1 <-
  interpola_cores_filtered_corit_PC1 %>%
  
  arrange(Age) %>%
  
  group_by(Age) %>%
  
  mutate(
    median_PC1 = median(value)
  ) %>%
  
  ungroup() %>%
  
  select(
    Age,
    median_PC1
  ) %>%
  
  distinct()

# ------------------------------------------------------------

# Regression and Age = 0 prediction

# ------------------------------------------------------------

median_values_PC1 <- median_value_PC1 %>%
  mutate(
    Age = Age / 1000
  )

lm_PC1 <- lm(
  median_PC1 ~ Age,
  data = median_values_PC1
)

predicted_PC1 <- predict(
  lm_PC1,
  newdata = data.frame(
    Age = 0
  )
)

predicted_row_PC1 <- data.frame(
  Age = 0,
  median_PC1 = predicted_PC1
)

median_values_PC1 <- rbind(
  predicted_row_PC1,
  median_values_PC1
) %>%
  arrange(Age)

# ------------------------------------------------------------

# PC1 plot

# ------------------------------------------------------------

PC1_plot <- ggplot(
  median_values_PC1,
  aes(
    x = Age,
    y = median_PC1
  )
) +
  
  geom_line(
    color = "#80B1D3",
    size = 1.2
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    size = 0.6
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      28,
      2
    ),
    limits = c(
      0,
      28
    )
  ) +
  
  scale_y_continuous(
    name = "PC1",
    position = "right",
    limits = c(
      -0.5,
      0.5
    )
  ) +
  
  theme_minimal() +
  
  theme(
    
    
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    panel.background = element_rect(
      fill = "white",
      color = "white"
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = "white"
    ),
    
    axis.text = element_text(
      size = 8
    ),
    
    axis.title = element_text(
      size = 8
    ),
    
    panel.border = element_rect(
      fill = NA,
      colour = "grey",
      size = 0.8
    ),
    
    axis.line = element_line(
      color = "grey",
      size = 0.8
    ),
    
    axis.ticks.length = unit(
      2,
      "mm"
    ),
    
    axis.ticks = element_line(
      color = "black",
      size = 1
    )
    
    
  )

# ============================================================

# 9. Combine all panels

# ============================================================

combined_plot <- (
  
  DNA_plot +
    plot_spacer() +
    
    Biomass_plot +
    plot_spacer() +
    
    TOC_plot +
    plot_spacer() +
    
    BR_Biomass_plot +
    plot_spacer() +
    
    ter_aqu_biomass +
    plot_spacer() +
    
    TJul_plot +
    plot_spacer() +
    
    PC1_plot
  
) +
  
  plot_layout(
    ncol = 1,
    
    
    heights = c(
      2,
      0.1,
      2,
      0.1,
      2,
      0.1,
      2,
      0.1,
      2,
      0.1,
      2,
      0.1,
      2
    ),
    
    guides = "collect"
    
    
  ) +
  
  plot_annotation(
    tag_levels = "a",
    
    
    theme = theme(
      plot.tag = element_text(
        size = 30,
        face = "bold"
      )
    )
    
    
  )

# ============================================================

# 10. Save Fig. 3

# ============================================================

ggsave(
  file.path(
    output_dir,
    "Fig3.png"
  ),
  combined_plot,
  width = 15,
  height = 20,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "Fig3.pdf"
  ),
  combined_plot,
  width = 15,
  height = 20
)

# ============================================================

# End of script

# ============================================================
