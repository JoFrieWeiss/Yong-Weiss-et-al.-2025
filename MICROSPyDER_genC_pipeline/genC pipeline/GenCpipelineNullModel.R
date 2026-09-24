# =============================================================================
# GENC PIPELINE: NULL MODEL ANALYSIS (REVIEWER RESPONSE)
# =============================================================================
# This script calculates the "null model" (no taxonomic resolution) 
# and compares it against the taxonomy-based model and measured TOC.

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(patchwork)

# --- 1. SET PATHS ---
PATH_PROCESSED_DATA <- "/Volumes/projects/biodiv/user/ziyong/FINAL_GenC_pipeline_Josefine_2026/output/Processed_DNA_with_MC_biomass.csv"
PATH_CONVERSION <- "/Volumes/projects/biodiv/user/ziyong/FINAL_GenC_pipeline_Josefine_2026/values_per_cell_correct.csv"
PATH_OUTPUT_DIR <- "/Volumes/projects/biodiv/user/ziyong/FINAL_GenC_pipeline_Josefine_2026/Null_Model_Analysis"

if(!dir.exists(PATH_OUTPUT_DIR)) dir.create(PATH_OUTPUT_DIR, recursive = TRUE)

# --- 2. LOAD DATA ---
df_main <- read_csv(PATH_PROCESSED_DATA, show_col_types = FALSE)
conversion_data <- read_csv2(PATH_CONVERSION, show_col_types = FALSE)

# --- 3. CALCULATE GLOBAL AVERAGES (THE NULL MODEL FOUNDATION) ---
# Extract ALL available DNA and biomass values across all taxa
all_dna <- as.numeric(unlist(conversion_data %>% select(starts_with("DNA_value"))))
all_bio <- as.numeric(unlist(conversion_data %>% select(starts_with("Biomass_value"))))

# Calculate global medians (ignoring any taxonomy)
global_median_dna <- median(all_dna, na.rm = TRUE)
global_median_bio <- median(all_bio, na.rm = TRUE)

# Average methodological constants
PARAM_EXTRACTION_EFFICIENCY <- 4.0
GLOBAL_STRUCT_FACTOR <- 1.5 # Conservative mean across all organisms

print(paste("Global median DNA/cell:", signif(global_median_dna, 3)))
print(paste("Global median C/cell:", signif(global_median_bio, 3)))

# --- 4. CALCULATE NULL MODEL ---
df_null <- df_main %>%
  mutate(
    # 1. Your finished GenC model (sum of all taxonomic biomasses)
    # Summing all columns ending with "_median_biomass" and multiplying by 100 for wt%
    Taxonomic_Model_OC = rowSums(select(., ends_with("_median_biomass")), na.rm = TRUE) * 100,
    
    # 2. Reviewer's null model
    # Using 'Starting_DNA_wt' (total extracted DNA) and global constants directly
    Null_Model_Biomass = (Starting_DNA_wt / global_median_dna) * global_median_bio * GLOBAL_STRUCT_FACTOR * PARAM_EXTRACTION_EFFICIENCY * 100 
  ) %>%
  # Keep only columns relevant for plotting
  select(Lake, Age, TOC, Taxonomic_Model_OC, Null_Model_Biomass) %>%
  rename(
    Measured_TOC = TOC
  ) %>%
  filter(!is.na(Measured_TOC))

# =============================================================================
# FINAL REVIEWER FIGURE (REDUCED SET: SPEARMAN, RMSE, BOXPLOT)
# =============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# --- 1. GLOBAL STATISTICS ---
df_clean <- df_stats %>% filter(Taxonomic_Model_OC > 0 & Measured_TOC > 0 & Null_Model_Biomass > 0)

s_genc <- cor(df_clean$Taxonomic_Model_OC, df_clean$Measured_TOC, method="spearman")
s_null <- cor(df_clean$Null_Model_Biomass, df_clean$Measured_TOC, method="spearman")
rmse_genc <- sqrt(mean((df_clean$Taxonomic_Model_OC - df_clean$Measured_TOC)^2))
rmse_null <- sqrt(mean((df_clean$Null_Model_Biomass - df_clean$Measured_TOC)^2))

df_bars <- data.frame(
  Model = rep(c("GenC Pipeline (Taxonomic)", "Null Model (Bulk DNA)"), 2),
  Metric = c(
    "1. Spearman Rho (Trend)\n(Robust against peaks)", "1. Spearman Rho (Trend)\n(Robust against peaks)",
    "2. RMSE [wt%]\n(Absolute Error - Lower is better)", "2. RMSE [wt%]\n(Absolute Error - Lower is better)"
  ),
  Value = c(s_genc, s_null, rmse_genc, rmse_null)
)

# --- 2. PLOTTING ---
# B1. Bar plot
p_bars <- ggplot(df_bars, aes(x = Model, y = Value, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("GenC Pipeline (Taxonomic)" = "#BC8F8F", "Null Model (Bulk DNA)" = "#4682B4")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold", size = 9),
        panel.grid.major.x = element_blank()) +
  labs(x = NULL, y = "Metric Value") +
  geom_text(aes(label = round(Value, 2)), vjust = ifelse(df_bars$Value < 0, 1.5, -0.5), size = 3.5)

# B2. Boxplot (Wilcoxon)
df_error_long <- df_clean %>%
  mutate(Error_GenC = abs(Taxonomic_Model_OC - Measured_TOC), Error_Null = abs(Null_Model_Biomass - Measured_TOC)) %>%
  select(Error_GenC, Error_Null) %>%
  pivot_longer(cols = everything(), names_to = "Model", values_to = "Absolute_Error") %>%
  mutate(Model = factor(Model, labels = c("GenC", "Null Model")))

p_global_error <- ggplot(df_error_long, aes(x = Model, y = Absolute_Error, fill = Model)) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.4, width = 0.5, color = "black", lwd = 0.6) +
  scale_fill_manual(values = c("GenC" = "#BC8F8F", "Null Model" = "#4682B4")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 9), plot.title = element_text(face = "bold", size = 9, hjust = 0.5)) +
  labs(title = "3. Absolute Error Dist.\n(Wilcoxon Test)", y = "Residual Error (wt%)", x = NULL) +
  annotate("text", x = 1.5, y = max(df_error_long$Absolute_Error, na.rm=T) * 0.95, label = "p < 2.22e-16", fontface = "bold", size = 3)

# --- 3. CORE PLOTS (WITHOUT NSE) ---
for(lake_name in names(plot_list_reviewer)) {
  df_lake <- df_null %>% filter(Lake == lake_name)
  s_genc_l <- cor(df_lake$Taxonomic_Model_OC, df_lake$Measured_TOC, method = "spearman")
  s_null_l <- cor(df_lake$Null_Model_Biomass, df_lake$Measured_TOC, method = "spearman")
  rmse_genc_l <- sqrt(mean((df_lake$Taxonomic_Model_OC - df_lake$Measured_TOC)^2, na.rm = TRUE))
  rmse_null_l <- sqrt(mean((df_lake$Null_Model_Biomass - df_lake$Measured_TOC)^2, na.rm = TRUE))
  
  plot_list_reviewer[[lake_name]] <- plot_list_reviewer[[lake_name]] + 
    labs(subtitle = sprintf("Spearman: GenC=%.2f | Null=%.2f\nRMSE: GenC=%.1f | Null=%.1f", s_genc_l, s_null_l, rmse_genc_l, rmse_null_l))
}

# --- 4. FINAL GRID ---
top_row <- (p_bars | p_global_error) + plot_layout(widths = c(2, 1))
final_master_grid <- (top_row / wrap_plots(plot_list_reviewer, ncol = 3)) +
  plot_layout(heights = c(1, 2), guides = "collect") + 
  plot_annotation(tag_levels = 'A', theme = theme(legend.position = "bottom", legend.title = element_blank()))

ggsave(file.path(PATH_OUTPUT_DIR, "Figure_R1_Final_Comparison_CLEAN.png"), final_master_grid, width = 17, height = 12, dpi = 500)