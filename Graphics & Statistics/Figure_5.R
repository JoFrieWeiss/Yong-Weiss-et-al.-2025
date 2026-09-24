
# =============================================================================
# Fig_5_plot
# Integrated comparison of aquatic and terrestrial biomass
# =============================================================================

# -----------------------------------------------------------------------------
# 0. LIBRARIES
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


# -----------------------------------------------------------------------------
# 1. PATHS
# -----------------------------------------------------------------------------

PATH_DATA_INPUT <- "data/Processed_DNA_with_MC_biomass.csv"
PATH_POLLEN <- "data/Pollendata/Algae pecentage csv.csv"
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
# 6. AQUATIC VS TERRESTRIAL BAR PLOT
# -----------------------------------------------------------------------------

Burial_rate_biomass_long <- Burial_rate_biomass %>%
  select(
    Startage,
    Biomass_Aquatic,
    Biomass_Terrestrial,
    core_name
  ) %>%
  pivot_longer(
    cols = starts_with("Biomass"),
    names_to = "Type",
    values_to = "Value"
  ) %>%
  mutate(
    Value = Value * 100
  ) %>%
  mutate(
    Type = factor(
      recode(
        Type,
        "Biomass_Aquatic" = "Aquatic",
        "Biomass_Terrestrial" = "Terrestrial"
      ),
      levels = c("Aquatic", "Terrestrial")
    )
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


box_biomass_ter_aqu_lake_bar <- Burial_rate_biomass_long %>%
  group_by(core_name, Type) %>%
  summarise(
    median_val = median(Value),
    q1 = quantile(Value, 0.25),
    q3 = quantile(Value, 0.75),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = core_name,
      y = median_val,
      fill = Type,
      color = Type
    )
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(ymin = q1, ymax = q3),
    position = position_dodge(width = 0.7),
    width = 0.2,
    size = 1.2,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Aquatic" = "#1f78b4",
      "Terrestrial" = "#33a02c"
    )
  ) +
  scale_color_manual(
    values = c(
      "Aquatic" = "#1f78b4",
      "Terrestrial" = "#33a02c"
    )
  ) +
  labs(
    y = "OC DNA-projected wt (%)",
    x = NULL,
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 7)) +
  base_theme +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box = "horizontal"
  )


# -----------------------------------------------------------------------------
# 7. LOAD POLLEN / BOTRYOCOCCUS DATA
# -----------------------------------------------------------------------------

clean_lake <- function(x) {
  x <- gsub("\u00A0", " ", x)
  x <- gsub("[^[:alnum:] ]", "", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}


Aquatic_data <- Burial_rate_biomass %>%
  select(
    Lake,
    Age,
    Biomass_Aquatic_percentage
  )

Aquatic_data$Lake_clean <- clean_lake(
  Aquatic_data$Lake
)


data <- read.csv(PATH_POLLEN)

data$Lake_clean <- clean_lake(data$Lake)


data <- data %>%
  mutate(
    Botryoccocus = ifelse(
      is.na(Botryoccocus),
      0,
      Botryoccocus
    ),
    total.pollen.sum = ifelse(
      is.na(total.pollen.sum),
      0,
      total.pollen.sum
    ),
    pct_Botry = ifelse(
      total.pollen.sum == 0,
      NA,
      Botryoccocus / total.pollen.sum * 100
    )
  )


# -----------------------------------------------------------------------------
# 8. PLOT EACH LAKE
# -----------------------------------------------------------------------------

lakes <- intersect(
  unique(Aquatic_data$Lake_clean),
  unique(data$Lake_clean)
)

lake_plots <- list()


for (lake in lakes) {
  
  aq <- Aquatic_data %>%
    filter(Lake_clean == lake) %>%
    arrange(Age)
  
  pol <- data %>%
    filter(Lake_clean == lake) %>%
    arrange(Age)
  
  if (nrow(aq) < 5 | nrow(pol) < 5) {
    next
  }
  
  
  # LOESS
  loess_model <- loess(
    pct_Botry ~ Age,
    data = pol,
    span = 0.3,
    na.action = na.exclude
  )
  
  aq$Botry_smooth <- predict(
    loess_model,
    newdata = data.frame(Age = aq$Age)
  )
  
  
  # Automatic scale
  scale_factor <- max(
    aq$Biomass_Aquatic_percentage,
    na.rm = TRUE
  ) /
    max(
      pol$pct_Botry,
      na.rm = TRUE
    )
  
  aq$Botry_smooth_scaled <-
    aq$Botry_smooth * scale_factor
  
  pol$pct_Botry_scaled <-
    pol$pct_Botry * scale_factor
  
  
  # Correlation
  smooth_clean <- aq %>%
    select(
      Biomass_Aquatic_percentage,
      Botry_smooth
    ) %>%
    drop_na()
  
  
  if (nrow(smooth_clean) > 3) {
    
    cor_smooth <- cor.test(
      smooth_clean$Biomass_Aquatic_percentage,
      smooth_clean$Botry_smooth,
      method = "spearman",
      exact = FALSE
    )
    
  } else {
    
    cor_smooth <- list(
      estimate = NA,
      p.value = NA
    )
    
  }
  
  
  get_sig_label <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("ns")
  }
  
  
  sig <- get_sig_label(
    cor_smooth$p.value
  )
  
  
  cor_label <- paste0(
    "r = ",
    round(cor_smooth$estimate, 2),
    "\n",
    "p = ",
    signif(cor_smooth$p.value, 3),
    " ",
    sig
  )
  
  
  # Plot
  p <- ggplot() +
    
    geom_line(
      data = aq,
      aes(
        Age,
        Biomass_Aquatic_percentage,
        color = "Biomass Aquatic"
      ),
      linewidth = 1.2
    ) +
    
    geom_point(
      data = pol,
      aes(
        Age,
        pct_Botry_scaled,
        color = "Botry (raw)"
      ),
      size = 1.2,
      alpha = 0.6
    ) +
    
    geom_line(
      data = aq,
      aes(
        Age,
        Botry_smooth_scaled,
        color = "Botry (LOESS)"
      ),
      linewidth = 1.5
    ) +
    
    scale_y_continuous(
      name = "Biomass Aquatic (%)",
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Botryococcus (%)"
      )
    ) +
    
    scale_color_manual(
      values = c(
        "Biomass Aquatic" = "steelblue",
        "Botry (raw)" = "darkblue",
        "Botry (LOESS)" = "darkblue"
      )
    ) +
    
    annotate(
      "text",
      x = max(aq$Age, na.rm = TRUE) * 0.65,
      y = max(
        aq$Biomass_Aquatic_percentage,
        na.rm = TRUE
      ) * 0.9,
      label = cor_label,
      size = 4
    ) +
    
    labs(
      title = lake,
      x = "Age"
    ) +
    
    theme_classic() +
    
    theme(
      legend.position = "top",
      legend.text = element_text(size = 10)
    )
  
  
  lake_plots[[lake]] <- p
}


# -----------------------------------------------------------------------------
# 9. COMBINE FIGURE 4.2
# -----------------------------------------------------------------------------

final_lakes_plot <- wrap_plots(
  lake_plots,
  ncol = 3
)


final_plot <- box_biomass_ter_aqu_lake_bar / final_lakes_plot +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Integrated Comparison of aquatic and terrestrial",
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
  file.path(PATH_OUTPUT_DIR, "Fig_5_plot.png"),
  final_plot,
  width = 36,
  height = 30,
  units = "cm",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(PATH_OUTPUT_DIR, "Fig_5_plot.pdf"),
  final_plot,
  width = 36,
  height = 30,
  units = "cm",
  device = cairo_pdf
)

