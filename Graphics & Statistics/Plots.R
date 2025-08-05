# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R SCRIPT FOR VISUALIZING LAKE SEDIMENT CORE DATA
#
# Description:
# This script loads and processes sediment core data to analyze biomass,
# DNA content, and other environmental proxies. It generates a series of
# figures (maps, time-series plots, PCA, RDA) to visualize the results.
#
# Author: Zijuan Yong
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## -----------------------------------------------------------------------------
## SECTION 1: LOAD LIBRARIES
## -----------------------------------------------------------------------------
# It's best practice to load all required libraries at the start of the script.

# Data manipulation and plotting
library(tidyverse) # Includes dplyr, ggplot2, tidyr, readr etc.
library(reshape2)  # For data reshaping (though pivot_longer/wider is preferred now)
library(cowplot)   # For combining ggplots
library(patchwork) # Alternative for combining ggplots
library(scales)    # For plot scaling functions (e.g., alpha, rescale)

# Statistical analysis and modeling
library(vegan)     # For ecological analyses like RDA
library(lmerTest)  # For linear mixed-effects models
library(effects)   # To visualize model effects
library(mgcv)      # For GAMs
library(corit)     # For interpolation (requires devtools)
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("EarthSystemDiagnostics/corit")

# Plotting enhancements
library(RColorBrewer) # For color palettes
library(ggrepel)      # To prevent text labels from overlapping in plots
library(ggtext)       # To render rich text (e.g., colored labels)
library(png)          # To read PNG files for plot composition
library(grid)         # For plot composition


## -----------------------------------------------------------------------------
## SECTION 2: CONFIGURATION AND DATA LOADING
## -----------------------------------------------------------------------------
# Define file paths.
# NOTE: Using relative paths (e.g., "data/my_file.csv") is recommended
# for better portability and use with RStudio Projects.
data_path <- "D:/document/DOC/AWI/Lake organic/data/Calculate_all_information_result_unclassfied.csv"
figure_path <- "D:/document/DOC/AWI/Lake organic/Figure/"
viridiplantae_path <- "D:/document/DOC/AWI/Lake organic/data/Viridiplantae_PCA_combin_all.csv"
rda_scores_path <- "D:/document/DOC/AWI/Lake organic/data/Plant_PCA_scores.csv"


# Create the figure directory if it doesn't exist
if (!dir.exists(figure_path)) {
  dir.create(figure_path, recursive = TRUE)
}

# --- Load and Prepare Main Dataset ---
burial_data <- read_csv(data_path)

# Replace NA values with 0
burial_data[is.na(burial_data)] <- 0

# Convert 'Age' to numeric in case it was loaded as a factor
burial_data$Age <- as.numeric(as.character(burial_data$Age))

# Standardize core names for consistency across all plots
burial_data <- burial_data %>%
  mutate(core_name = case_when(
    grepl("Ilirney", Core, ignore.case = TRUE) ~ "Ilirney",
    grepl("Lele", Core, ignore.case = TRUE) ~ "Lele",
    grepl("Lama", Core, ignore.case = TRUE) ~ "Lama",
    grepl("Btoko", Core, ignore.case = TRUE) ~ "BToko",
    grepl("Salmon", Core, ignore.case = TRUE) ~ "Salmon",
    grepl("Ulu", Core, ignore.case = TRUE) ~ "Ulu",
    TRUE ~ as.character(Core) # Keep original name if no match
  ))

# Define a function to assign time slices based on age
assign_timeslice <- function(age, interval = 2, max_age = 60) {
  timeslices <- seq(0, max_age, by = interval)
  slice_index <- findInterval(age, timeslices, rightmost.closed = TRUE)
  start <- timeslices[slice_index]
  end <- timeslices[slice_index] + interval
  if (age >= max_age) return(paste0(max_age, "+"))
  return(paste0(start, "-", end))
}

# Apply the timeslice function to create time bins
burial_data <- burial_data %>%
  mutate(
    timeslice = sapply(Age, assign_timeslice, interval = 2, max_age = 60),
    Startage = as.integer(sub("-.*", "", timeslice))
  )

# Define shared plotting aesthetics
core_colors <- c(
  "Ilirney" = "brown", "Lele" = "forestgreen", "Lama" = "royalblue4",
  "BToko" = "#87CEFF", "Salmon" = "darkolivegreen", "Ulu" = "orange"
)
core_levels <- c("BToko", "Ulu", "Lama", "Ilirney", "Salmon", "Lele")

base_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey50", linewidth = 0.7),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

## -----------------------------------------------------------------------------
## SECTION 3: EXPLORATORY DATA ANALYSIS (EDA)
## -----------------------------------------------------------------------------
# ... (This section is the same as in the previous response) ...


## -----------------------------------------------------------------------------
## FIGURE 2: Bar plot of OC burial rate per lake core
## -----------------------------------------------------------------------------
# ... (This section is the same as in the previous response) ...


## -----------------------------------------------------------------------------
## FIGURE 3: Composite time-series plots
## -----------------------------------------------------------------------------
# ... (This section is the same as in the previous response) ...


## -----------------------------------------------------------------------------
## FIGURE 4: Comparison of TOC and DNA-projected OC
## -----------------------------------------------------------------------------
# ... (This section is the same as in the previous response) ...


## -----------------------------------------------------------------------------
## FIGURE 5: DNA vs. Biomass Weight Percentage by Taxon
## -----------------------------------------------------------------------------
# This figure compares the contribution of different taxonomic groups based on
# their raw DNA weight vs. their projected biomass weight. It uses a
# pseudo-logarithmic scale to accommodate a wide range of values.

# Define taxon order and colors for consistent plotting
taxon_order <- c(
  "Aquatic algae", "Archaea aquatic", "Bacteria aquatic", "Fungi aquatic",
  "Metazoa aquatic", "Viridiplantae aquatic", "Archaea terrestrial",
  "Bacteria terrestrial", "Fungi terrestrial", "Metazoa terrestrial",
  "Viridiplantae woody", "Viridiplantae non woody", "Viruses"
)
color_values <- c(
  "#ADD8E6", "#00BFFF", "#4169E1", "#0000FF", "#00008B", "royalblue4",
  "#98FB98", "#50EE50", "#32CD32", "#228B22", "#336433", "#006400", "#003300"
)

# --- Panel A: DNA vs. Biomass (Absolute Weight %) ---

# Prepare DNA data
df_dna <- burial_data %>%
  select(starts_with("Aquatic_algae_percentage_Start_wt"):starts_with("Viruses_percentage_Start_wt")) %>%
  rename_with(~ str_remove(., "_percentage_Start_wt")) %>%
  rename_with(~ str_replace_all(., "_", " ")) %>%
  pivot_longer(everything(), names_to = "Taxon", values_to = "Value") %>%
  mutate(Value = Value * 100, Source = "DNA weight percentage")

# Prepare Biomass data
df_biomass <- burial_data %>%
  select(starts_with("median_biomass_Aquatic_algae"):starts_with("median_biomass_Viruses")) %>%
  rename_with(~ str_remove(., "median_biomass_")) %>%
  rename_with(~ str_replace_all(., "_", " ")) %>%
  pivot_longer(everything(), names_to = "Taxon", values_to = "Value") %>%
  mutate(Value = Value * 100, Source = "Biomass weight percentage")

# Combine and summarize data
df_combined_summary <- bind_rows(df_dna, df_biomass) %>%
  group_by(Taxon, Source) %>%
  summarise(
    median = median(Value, na.rm = TRUE),
    q1 = quantile(Value, 0.25, na.rm = TRUE),
    q3 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order))

# A pseudo-log transformation is used to visualize data spanning several orders of magnitude,
# including very small values, which would be compressed on a standard log scale.
pseudo_log_trans <- function(sigma = 1e-6) {
  trans_new("pseudo_log",
            transform = function(x) asinh(x / (2 * sigma)),
            inverse = function(x) 2 * sigma * sinh(x))
}

plot_a_fig5 <- ggplot(df_combined_summary, aes(x = Taxon, y = median, group = Source)) +
  # Use ifelse to create hollow bars (DNA) and filled bars (Biomass)
  geom_col(
    aes(fill = ifelse(Source == "Biomass weight percentage", as.character(Taxon), NA),
        color = as.character(Taxon),
        linewidth = ifelse(Source == "DNA weight percentage", 1.2, 0.5)),
    position = position_dodge(width = 0.7), width = 0.6
  ) +
  geom_errorbar(aes(ymin = q1, ymax = q3, color = Taxon),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 1.2) +
  scale_fill_manual(values = setNames(alpha(color_values, 0.6), taxon_order), na.value = "transparent") +
  scale_color_manual(values = setNames(color_values, taxon_order), guide = "none") +
  scale_linewidth_identity() +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-6),
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(y = "Weight percentage (wt%, log10 scale)", x = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# --- Panel B: Biomass Composition (Relative %) ---

# Prepare data for relative biomass contribution
df_relative_biomass <- burial_data %>%
  select(matches("median_biomass_.*_percent$")) %>%
  rename_with(~ str_remove(., "median_biomass_")) %>%
  rename_with(~ str_remove(., "_percent")) %>%
  rename_with(~ str_replace_all(., "_", " ")) %>%
  pivot_longer(everything(), names_to = "Taxon", values_to = "Percentage") %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order))

# Summarize data
df_relative_summary <- df_relative_biomass %>%
  group_by(Taxon) %>%
  summarise(
    median = median(Percentage, na.rm = TRUE),
    q1 = quantile(Percentage, 0.25, na.rm = TRUE),
    q3 = quantile(Percentage, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

plot_b_fig5 <- ggplot(df_relative_summary, aes(x = Taxon, y = median, fill = Taxon, color = Taxon)) +
  geom_col(width = 0.6, linewidth = 1.2) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.2, linewidth = 1.2) +
  scale_x_discrete(limits = taxon_order) +
  scale_y_continuous(name = "Contribution to total OC DNA-projected wt (%)") +
  scale_fill_manual(values = alpha(color_values, 0.6)) +
  scale_color_manual(values = color_values) +
  coord_cartesian(ylim = c(0, 42)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

# --- Combine for Figure 5 ---
combined_plot_fig5 <- plot_a_fig5 / plot_b_fig5 +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(figure_path, "fig5_DNA_vs_Biomass.png"), combined_plot_fig5, width = 12, height = 10)
ggsave(file.path(figure_path, "fig5_DNA_vs_Biomass.pdf"), combined_plot_fig5, width = 12, height = 10)


## -----------------------------------------------------------------------------
## FIGURE 6 & S11: Time-series of Biomass Contribution by Taxon
## -----------------------------------------------------------------------------
# These plots show the change in biomass contribution for each taxonomic group
# over time, with each group in its own panel.
# Fig 6 uses relative percentages. Fig S7 uses absolute weight percentages.

# Generic function to create the faceted time-series plot
create_faceted_timeseries <- function(data, value_col, y_label) {
  # Prepare data
  data_long <- data %>%
    select(Startage, matches(value_col)) %>%
    pivot_longer(
      cols = -Startage,
      names_to = "Type",
      values_to = "Value"
    ) %>%
    mutate(
      Type = str_remove_all(Type, "median_biomass_|_percent"),
      Type = str_replace_all(Type, "_", " "),
      Type = factor(Type, levels = taxon_order)
    )
  
  # Summarize data by time slice
  summary_data <- data_long %>%
    group_by(Startage, Type) %>%
    summarise(
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Startage = factor(Startage))
  
  # Create plot
  ggplot(summary_data, aes(x = Startage, y = Median, fill = Type, color = Type)) +
    geom_bar(stat = "identity", width = 0.7, position = "dodge") +
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, linewidth = 1, position = "dodge") +
    facet_wrap(~Type, scales = "free_y", ncol = 1, strip.position = "top") +
    scale_fill_manual(values = alpha(color_values, 0.6), breaks = taxon_order) +
    scale_color_manual(values = color_values, breaks = taxon_order) +
    labs(x = "Age (kyr)", y = y_label) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 10, hjust = 0),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# --- Generate Figure 6 (Relative Percentages) ---
plot_fig6 <- create_faceted_timeseries(
  data = burial_data %>% filter(Age < 30),
  value_col = "_percent$",
  y_label = "Contribution to total OC DNA-projected wt (%)"
)
ggsave(file.path(figure_path, "fig6_percent_wt.png"), plot_fig6, width = 8, height = 20)
ggsave(file.path(figure_path, "fig6_percent_wt.pdf"), plot_fig6, width = 8, height = 20)

# --- Generate Figure S7 (Absolute Weight Percentages) ---
plot_figs7 <- create_faceted_timeseries(
  data = burial_data %>% filter(Age < 30),
  value_col = "median_biomass_((?!percent).)*$", # Regex for columns without 'percent'
  y_label = "OC DNA-projected weight (wt%)"
)
ggsave(file.path(figure_path, "FigS7_wt_per_group.png"), plot_figs7, width = 8, height = 20)
ggsave(file.path(figure_path, "FigS7_wt_per_group.pdf"), plot_figs7, width = 8, height = 20)


## -----------------------------------------------------------------------------
## FIGURE S5: TOC vs DNA-projected OC Time-series per Lake
## -----------------------------------------------------------------------------
# This plot shows trends of TOC and DNA-projected OC over the last 60 kyr,
# faceted for each individual lake core.

# Prepare data
data_long_s3 <- burial_data %>%
  mutate(Total_median_biomass = Total_median_biomass * 100) %>%
  select(Age, core_name, TOC, Total_median_biomass) %>%
  pivot_longer(
    cols = c(TOC, Total_median_biomass),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(Type = recode(Type, "Total_median_biomass" = "OC DNA-projected"))

# Create plot
plot_s3 <- ggplot(data_long_s3, aes(x = Age, y = Value, color = Type)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2) +
  facet_wrap(~core_name, scales = "free_y") +
  scale_x_continuous(name = "Age (kyr)", limits = c(0, 60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(name = "Weight percentage (%)") +
  scale_color_manual(values = c("TOC" = "#5d2c04", "OC DNA-projected" = "#BC8F8F")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", legend.title = element_blank(), panel.grid = element_blank())

ggsave(file.path(figure_path, "Figs3_TOC_vs_DNA_per_lake.png"), plot_s3, width = 12, height = 8)
ggsave(file.path(figure_path, "Figs3_TOC_vs_DNA_per_lake.pdf"), plot_s3, width = 12, height = 8)


## -----------------------------------------------------------------------------
## FIGURE S6: LMM/GLM Results for Burial Rates
## -----------------------------------------------------------------------------
# This section fits Linear Mixed-Effects Models to determine the drivers of
# total, aquatic, and terrestrial OC burial rates.

# Prepare data for modeling
env_std <- burial_data %>%
  select(Number, Age, Lake, Sediment.Rate, PC1, PC2, Pann, TJul) %>%
  group_by(Lake) %>%
  mutate(across(c(Age, Sediment.Rate, PC1, PC2, Pann, TJul),
                ~as.vector(scale(.)), .names = "scaled_{.col}")) %>%
  ungroup()

# --- Generic Function to Create Model Plots ---
create_glm_plot <- function(response_var, data, env_data, plot_label, main_title) {
  
  # Prepare response data and combine with predictors
  response_df <- data %>%
    select(Number, Age, all_of(response_var)) %>%
    mutate(Response = log10(.data[[response_var]]))
  
  model_data <- inner_join(response_df, env_data, by = c("Number", "Age"))
  
  # Fit the model
  fit <- lmer(Response ~ scaled_Age + scaled_Sediment.Rate + scaled_PC1 +
                scaled_PC2 + scaled_Pann + scaled_TJul + (1 | Lake),
              data = model_data)
  
  # 1. Create diagnostic plot and save to a temporary file
  temp_png <- tempfile(fileext = ".png")
  png(temp_png, width = 800, height = 800, res = 150)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  plot(fit)
  dev.off()
  
  # 2. Create effect plots
  effect_plots <- plot(effect("scaled_Age", fit), plot = FALSE, main=NULL) +
    plot(effect("scaled_Sediment.Rate", fit), plot = FALSE, main=NULL) +
    plot(effect("scaled_PC1", fit), plot = FALSE, main=NULL) +
    plot(effect("scaled_PC2", fit), plot = FALSE, main=NULL) +
    plot(effect("scaled_Pann", fit), plot = FALSE, main=NULL) +
    plot(effect("scaled_TJul", fit), plot = FALSE, main=NULL) +
    plot_layout(ncol = 3)
  
  # 3. Combine diagnostic and effect plots
  residual_grob <- rasterGrob(readPNG(temp_png), interpolate = TRUE)
  
  combined <- plot_grid(residual_grob, effect_plots, ncol = 2,
                        labels = plot_label, label_size = 18,
                        rel_widths = c(1, 1.5))
  
  # Add a main title
  final_plot <- ggdraw() +
    draw_plot(combined) +
    draw_label(main_title, fontface = 'bold', size = 16, x = 0.5, y = 0.98)
  
  return(final_plot)
}

# --- Generate plots for each response variable ---
plot_total_br <- create_glm_plot("BR_biomass", burial_data, env_std, c("A", "B"), "Total OC DNA-projected Burial Rate")
plot_aquatic_br <- create_glm_plot("BR_biomass_Aquatic", burial_data, env_std, c("C", "D"), "Aquatic OC DNA-projected Burial Rate")
plot_terrestrial_br <- create_glm_plot("BR_biomass_terrestrial", burial_data, env_std, c("E", "F"), "Terrestrial OC DNA-projected Burial Rate")

# --- Combine all three model plots ---
final_s6_plot <- plot_grid(plot_total_br, plot_aquatic_br, plot_terrestrial_br, ncol = 1)

ggsave(file.path(figure_path, "Fig_S6_GLM_BR.png"), final_s8_plot, width = 16, height = 24, dpi = 300)
ggsave(file.path(figure_path, "Fig_S6_GLM_BR.pdf"), final_s8_plot, width = 16, height = 24) 

            
## -----------------------------------------------------------------------------
## FIGURE S7 & S8: PCA of Viridiplantae Composition
## -----------------------------------------------------------------------------

# --- Load and Prepare Viridiplantae Data ---
virid_data <- read_csv(viridiplantae_path) %>%
  mutate(
    Lake_group = factor(case_when(
      grepl("Ilirney", Lake, ignore.case = TRUE) ~ "Ilirney",
      grepl("Lele", Lake, ignore.case = TRUE) ~ "Lele",
      grepl("Lama", Lake, ignore.case = TRUE) ~ "Lama",
      grepl("Btoko", Lake, ignore.case = TRUE) ~ "BToko",
      grepl("Ulu", Lake, ignore.case = TRUE) ~ "Ulu",
      grepl("Salmon", Lake, ignore.case = TRUE) ~ "Salmon"
    ), levels = core_levels)
  )

# --- Perform PCA (using RDA with no constraints) ---
# Select species data columns
species_data <- virid_data %>% select(3:30)

# Hellinger transformation is often used for species abundance data to down-weight rare species
species_hellinger <- decostand(species_data, method = "hellinger")
species_hellinger[is.na(species_hellinger)] <- 0

# Perform PCA/RDA
rda_result <- rda(species_hellinger)

# Save the scores for later use if needed
site_scores_raw <- scores(rda_result, display = "sites")
# write.csv(site_scores_raw, rda_scores_path, row.names = FALSE)


# --- FIGURE S7: PCA Biplot ---

# Extract and normalize scores for plotting
site_scores <- scores(rda_result, display = "sites")
species_scores <- scores(rda_result, display = "species")
explained_variance <- summary(rda_result)$cont$importance[2, 1:2] * 100

# Function to normalize scores to a [-1, 1] range for visualization
normalize_range <- function(x) (x - min(x)) / (max(x) - min(x)) * 2 - 1

# Prepare data for plotting
site_scores_df <- as.data.frame(apply(site_scores, 2, normalize_range)) %>%
  mutate(
    Lake_group = virid_data$Lake_group,
    Age = virid_data$Age,
    Startage = as.integer(sub("-.*", "", sapply(Age, assign_timeslice)))
  )

species_scores_df <- as.data.frame(apply(species_scores, 2, normalize_range)) %>%
  mutate(species = rownames(.))

# Calculate mean position for each time slice per lake
mean_scores <- site_scores_df %>%
  group_by(Lake_group, Startage) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = 'drop')

# Create colored labels for the legend using ggtext
legend_labels <- paste0("<span style='color:", core_colors, ";'>", names(core_colors), "</span>")

plot_s7 <- ggplot(mean_scores, aes(x = PC1, y = PC2, color = Lake_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_text_repel(aes(label = Startage), size = 4, max.overlaps = Inf) +
  geom_text_repel(data = species_scores_df, aes(x = PC1, y = PC2, label = species),
                  inherit.aes = FALSE, size = 4, max.overlaps = Inf) +
  scale_color_manual(values = core_colors, labels = legend_labels) +
  labs(
    title = "PCA of Viridiplantae Composition",
    x = sprintf("PC1 (%.1f%%)", explained_variance[1]),
    y = sprintf("PC2 (%.1f%%)", explained_variance[2]),
    color = "Lake"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.text = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path(figure_path, "FigS7_Viridiplantae_PCA.pdf"), plot_s4, width = 10, height = 7)
ggsave(file.path(figure_path, "FigS7_Viridiplantae_PCA.png"), plot_s4, width = 10, height = 7)


# --- FIGURE S8: PC1 Scores vs. Time ---
plot_s8 <- ggplot(site_scores_df, aes(x = Age, y = PC1)) +
  geom_line(color = "#80B1D3", linewidth = 1) +
  geom_point(color = "#80B1D3", size = 2, alpha = 0.7) +
  facet_wrap(~Lake_group, ncol = 2) +
  scale_x_continuous(name = "Age (kyr)", breaks = seq(0, 60, 10)) +
  scale_y_continuous(name = "PC1 Score", limits = c(-1, 1)) +
  labs(title = "Time-series of Viridiplantae PC1 Scores") +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_path, "FigS8_PCA_scores.png"), plot_s5, width = 10, height = 8)
ggsave(file.path(figure_path, "FigS8_PCA_scores.pdf"), plot_s5, width = 10, height = 8)


## -----------------------------------------------------------------------------
## FIGURE S9: Time-series of Salicaceae and Rosaceae
## -----------------------------------------------------------------------------
# This plot shows the percentage of two key plant families over the last 30 kyr.

virid_data_30k <- read_csv(viridiplantae_path) %>%
  filter(Age < 30) %>%
  mutate(Startage = as.integer(sub("-.*", "", sapply(Age, assign_timeslice))))

# Generic plotting function for a single family
plot_family_timeseries <- function(data, family_name, color, show_x_axis = TRUE) {
  summary_data <- data %>%
    group_by(Startage) %>%
    summarise(
      Median = median(.data[[family_name]], na.rm = TRUE),
      Q1 = quantile(.data[[family_name]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[family_name]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(summary_data, aes(x = as.factor(Startage), y = Median)) +
    geom_col(fill = alpha(color, 0.6), width = 0.7) +
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, color = color, linewidth = 1) +
    labs(y = paste(family_name, "(%)"), x = if (show_x_axis) "Age (kyr)" else NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  return(p)
}

# Create and combine plots
plot_salicaceae <- plot_family_timeseries(virid_data_30k, "Salicaceae", "#336633", show_x_axis = FALSE)
plot_rosaceae <- plot_family_timeseries(virid_data_30k, "Rosaceae", "#228B22", show_x_axis = TRUE)

combined_s9 <- plot_salicaceae / plot_rosaceae

ggsave(file.path(figure_path, "FigS9_Families.png"), combined_s11, width = 8, height = 6)
ggsave(file.path(figure_path, "FigS9_Families.pdf"), combined_s11, width = 8, height = 6)
            
## -----------------------------------------------------------------------------
## FIGURE S10: Redundancy Analysis (RDA)
## -----------------------------------------------------------------------------
# This RDA explores the relationship between taxonomic composition (all groups)
# and environmental variables.

# Prepare species data (using relative biomass percentages)
species_rda <- burial_data %>%
  select(matches("median_biomass_.*_percent$")) %>%
  rename_with(~ str_remove(., "median_biomass_|_percent")) %>%
  rename_with(~ str_replace_all(., "_", " "))
# Fourth-root transformation to stabilize variance and down-weight dominant taxa
species_rda_scaled <- sqrt(sqrt(species_rda))
species_rda_scaled[is.na(species_rda_scaled)] <- 0

# Prepare environmental data (standardize within each lake)
env_rda <- burial_data %>%
  select(Lake, Age, PC1, PC2, Pann, TJul, Sediment.Rate) %>%
  group_by(Lake) %>%
  mutate(across(-Lake, ~scale(.)[,1])) %>% # scale returns a matrix
  ungroup()

# Perform RDA
rda_env_result <- rda(species_rda_scaled ~ Age + Sediment.Rate + PC1 + PC2 + Pann + TJul, data = env_rda)
rda_summary <- summary(rda_env_result)
rda_explained <- rda_summary$cont$importance[2, 1:2] * 100

# Extract scores and normalize for plotting
site_scores_rda <- normalize_range(scores(rda_env_result, display = "sites")[,1:2])
species_scores_rda <- normalize_range(scores(rda_env_result, display = "species")[,1:2])
env_scores_rda <- scores(rda_env_result, display = "bp")

# Prepare data for ggplot
rda_sites_df <- as.data.frame(site_scores_rda) %>%
  mutate(
    Lake_group = factor(burial_data$core_name, levels = core_levels),
    Startage = burial_data$Startage
  )

rda_species_df <- as.data.frame(species_scores_rda) %>% mutate(label = rownames(.))
rda_env_df <- as.data.frame(env_scores_rda) %>% mutate(label = rownames(.))

# Get mean positions for time slice labels
rda_mean_scores <- rda_sites_df %>%
  group_by(Lake_group, Startage) %>%
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2), .groups = "drop")

plot_s10 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  # Environmental variable arrows
  geom_segment(data = rda_env_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "inches")), color = "black") +
  geom_text_repel(data = rda_env_df, aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = label), size = 5) +
  # Species labels
  geom_text_repel(data = rda_species_df, aes(x = RDA1, y = RDA2, label = label), size = 4) +
  # Time slice labels
  geom_text_repel(data = rda_mean_scores, aes(x = RDA1, y = RDA2, label = Startage, color = Lake_group), size = 4) +
  scale_color_manual(values = core_colors, labels = legend_labels) +
  labs(
    title = "Redundancy Analysis (RDA)",
    x = sprintf("RDA1 (%.1f%%)", rda_explained[1]),
    y = sprintf("RDA2 (%.1f%%)", rda_explained[2]),
    color = "Lake"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.text = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path(figure_path, "FigS10_RDA.pdf"), plot_s6, width = 12, height = 9)
ggsave(file.path(figure_path, "FigS10_RDA.png"), plot_s6, width = 12, height = 9)

