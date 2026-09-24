# =============================================================================
# GENC PIPELINE PART 3: BIOMASS RECONSTRUCTION
# =============================================================================
## please change path for windows in line 27, 28, 29
# ==============================================================================

# --- 0. LIBRARIES ---
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(stringr)
if(!require(sensitivity)) install.packages("sensitivity")
library(sensitivity)
if(!require(patchwork)) install.packages("patchwork")
library(patchwork)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)

# =============================================================================
# >>> USER CONFIGURATION <<<
# =============================================================================

PATH_DATA_INPUT <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/output/Processed_DNA_habitat_percentage_SPLIT_UNCLASSIFIED.csv"
PATH_CONVERSION <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/values_per_cell_correct.csv"
PATH_OUTPUT_DIR <- "D:/document/DOC/AWI/Lake organic/FINAL_GenC_pipeline_Josefine_2026/MC_Results_Final_Robust"

# --- B. SIMULATION SETTINGS ---
N_SIMS_MC <- 10000      # High N stabilizes the Median against outliers
N_SOBOL   <- 2000       

# --- C. CORRECTION FACTORS  ---
PARAM_EXTRACTION_EFFICIENCY <- 4.0 
PARAM_STRUCT_BACTERIA <- 2.5   
PARAM_STRUCT_ALGAE    <- 1.5   
PARAM_STRUCT_WOODY    <- 4.0  
PARAM_STRUCT_PLANT    <- 1.0   


# --- D. PERTURBATION RANGES (STRESS TEST) ---
# Small Steps: Test methodological uncertainty (+/- 10%, 25%)
PERT_STEPS_METHOD <- c(0.75, 0.90, 1.10, 1.25) 
# Large Steps: Test fundamental biological shifts (Half/Double)
PERT_STEPS_BIO    <- c(0.50, 2.00)             

# =============================================================================
# END OF CONFIGURATION
# =============================================================================

if(!dir.exists(PATH_OUTPUT_DIR)) dir.create(PATH_OUTPUT_DIR, recursive = TRUE)
simulation_points_cache <- list()

# --- PLOTTING COLORS ---

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
  "Viruses" )

color_values <- c(
  "#ADD8E6", "#00BFFF", "#4169E1", "#0000FF", "#00008B", "royalblue4",
  "#98FB98", "#50EE50", "#32CD32", "#228B22", "#336433", "#006400", "#003300"
)
names(color_values) <- taxon_order[1:length(color_values)]
color_values_alpha <- alpha(color_values, 0.6)

# -----------------------------------------------------------------------------
# 1. LOAD DATA (ROBUST VERSION)
# -----------------------------------------------------------------------------
print(">>> LOADING DATA...")

if(!file.exists(PATH_CONVERSION)) stop(paste("Conversion file not found:", PATH_CONVERSION))

if(file.exists(PATH_DATA_INPUT)) {
  DNA_global <- read_csv(PATH_DATA_INPUT, show_col_types = FALSE)
} else if (exists("combine_data_calc")) {
  print("File not found, using 'combine_data_calc' from environment.")
  DNA_global <- combine_data_calc
} else {
  stop("Input Data not found.")
}

# --- STEP A: FIX DNA WEIGHT NAMES ---
# Convert Part 2 names (_Start_wt) to Part 3 names (_DNA_weight)
names(DNA_global) <- names(DNA_global) %>% 
  gsub("_Start_wt", "_DNA_weight", .) %>% 
  gsub("_percentage", "", .)

# --- STEP B: DYNAMIC COLUMN REPAIR (Age and Core) ---
# Check if Age.x exists, if not use Age, if not use join_Age
if("Age.x" %in% names(DNA_global)) {
  DNA_global$Age <- DNA_global$Age.x
} else if ("join_Age" %in% names(DNA_global)) {
  DNA_global$Age <- DNA_global$join_Age
} 
# Note: If it's already called "Age", you don't need to do anything.

# Same for Core
if("Core.x" %in% names(DNA_global)) {
  DNA_global$Core <- DNA_global$Core.x
} else if ("join_Core" %in% names(DNA_global)) {
  DNA_global$Core <- DNA_global$join_Core
}

# --- STEP C: FINAL VERIFICATION ---
if(!"Age" %in% names(DNA_global)) {
  # Fallback: Find any column that contains "Age"
  existing_age_col <- names(DNA_global)[grepl("Age", names(DNA_global), ignore.case = TRUE)][1]
  if(!is.na(existing_age_col)) DNA_global$Age <- DNA_global[[existing_age_col]]
}



values_per_cell_global <- read_csv2(PATH_CONVERSION, show_col_types = FALSE)
core_list <- sort(unique(DNA_global$Number[!is.na(DNA_global$Number)]))

message("Data loaded and Age/Core columns verified.")

# =============================================================================
# PART 2: FUNCTIONS
# =============================================================================

# --- A. Monte Carlo Simulation  ---
# We stabilize the formula by using Beta-Distributions (Priors) that prevent
# unrealistic genome sizes for deep biosphere organisms.
run_monte_carlo <- function(data, dna_col_name, organism_name, conversion_data, n_sims, 
                            ext_factor = PARAM_EXTRACTION_EFFICIENCY, 
                            struct_factors = list(bac=PARAM_STRUCT_BACTERIA, alg=PARAM_STRUCT_ALGAE, woody=PARAM_STRUCT_WOODY, plant=PARAM_STRUCT_PLANT),
                            dna_scaling = 1.0, bio_scaling = 1.0) {
  
  factors <- conversion_data %>% filter(Organism == organism_name)
  if(nrow(factors) == 0 || !(dna_col_name %in% names(data))) return(NULL)
  
  dna_pool <- as.numeric(unlist(factors %>% select(starts_with("DNA_value")))); dna_pool <- dna_pool[!is.na(dna_pool)]
  bio_pool <- as.numeric(unlist(factors %>% select(starts_with("Biomass_value")))); bio_pool <- bio_pool[!is.na(bio_pool)]
  
  weights <- as.numeric(as.character(data[[dna_col_name]])); weights[is.na(weights)] <- 0
  if(length(weights) == 0) return(NULL)
  
  sim_results <- matrix(NA, nrow = length(weights), ncol = n_sims)
  d_min <- min(dna_pool); d_max <- max(dna_pool)
  b_min <- min(bio_pool); b_max <- max(bio_pool)
  CARBON_CONVERSION <- 1 
  
  # --- LOGIC BLOCKS WITH STABILIZING PRIORS ---
  if (grepl("Bacteria", organism_name, ignore.case = TRUE)) {
    val_min <- 3.00e-14; val_max <- 1.18e-12
    # Diagnostic
    diag_dna <- (d_min + rbeta(100, 1, 20) * (d_max - d_min)) * dna_scaling
    diag_bio <- ((val_min + rbeta(100, 1, 2) * (val_max - val_min)) * bio_scaling) * struct_factors$bac * ext_factor
    
    for (i in 1:n_sims) {
      # Prior: Beta(1,10) skews towards small genomes -> Stabilizes against overestimation
      rand_dna <- (d_min + rbeta(length(weights), 1, 20) * (d_max - d_min)) * dna_scaling
      c_content <- (val_min + rbeta(length(weights), 1, 2) * (val_max - val_min)) * bio_scaling
      sim_results[, i] <- (weights / rand_dna) * c_content * struct_factors$bac * ext_factor * CARBON_CONVERSION
    }
    
  } else if (grepl("Algae|Phytoplankton", organism_name, ignore.case = TRUE) || grepl("Aquatic_algae", dna_col_name)) {
    diag_dna <- (d_min + rbeta(100, 1, 10) * (d_max - d_min)) * dna_scaling
    diag_bio <- ((b_min + rbeta(100, 10, 1) * (b_max - b_min)) * bio_scaling) * struct_factors$alg * ext_factor
    for (i in 1:n_sims) {
      rand_dna <- (d_min + rbeta(length(weights), 1, 10) * (d_max - d_min)) * dna_scaling
      rand_bio <- (b_min + rbeta(length(weights), 10, 1) * (b_max - b_min)) * bio_scaling
      sim_results[, i] <- (weights / rand_dna) * rand_bio * struct_factors$alg * ext_factor * CARBON_CONVERSION
    }
    
  } else if (grepl("Woody plant", organism_name, ignore.case = TRUE)) {
    diag_dna <- (d_min + rbeta(100, 1, 2) * (d_max - d_min)) * dna_scaling
    diag_bio <- ((b_min + rbeta(100, 4, 1) * (b_max - b_min)) * bio_scaling) * struct_factors$woody * ext_factor
    for (i in 1:n_sims) {
      rand_dna <- (d_min + rbeta(length(weights), 1, 2) * (d_max - d_min)) * dna_scaling
      rand_bio <- (b_min + rbeta(length(weights), 4, 1) * (b_max - b_min)) * bio_scaling
      sim_results[, i] <- (weights / rand_dna) * rand_bio * struct_factors$woody * ext_factor * CARBON_CONVERSION
    }
    
  } else if (grepl("Plant|Viridiplantae", dna_col_name, ignore.case = TRUE)) {
    diag_dna <- (d_min + rbeta(100, 1, 1) * (d_max - d_min)) * dna_scaling
    diag_bio <- (runif(100, b_min, b_max) * bio_scaling) * struct_factors$plant * ext_factor
    for (i in 1:n_sims) {
      rand_dna <- (d_min + rbeta(length(weights), 1, 1) * (d_max - d_min)) * dna_scaling
      rand_bio <- (runif(length(weights), b_min, b_max) * bio_scaling)
      sim_results[, i] <- (weights / rand_dna) * rand_bio * struct_factors$plant * ext_factor * CARBON_CONVERSION
    }
    
  } else {
    diag_dna <- runif(100, d_min, d_max) * dna_scaling
    diag_bio <- runif(100, b_min, b_max) * bio_scaling * 1.0 * ext_factor
    for (i in 1:n_sims) {
      rand_dna <- runif(length(weights), d_min, d_max) * dna_scaling
      rand_bio <- runif(length(weights), b_min, b_max) * bio_scaling
      sim_results[, i] <- (weights / rand_dna) * rand_bio * 1.0 * ext_factor * CARBON_CONVERSION
    }
  }
  
  if(!is.null(diag_dna) && length(simulation_points_cache) < 2000) { 
    disp_name <- organism_name; if(grepl("Bacteria", organism_name)) disp_name <- "Bacteria"
    df_diag <- data.frame(Group = disp_name, DNA_Factor = diag_dna, Cell_Biomass = diag_bio)
    simulation_points_cache <<- bind_rows(simulation_points_cache, df_diag)
  }
  
  data.frame(
    Age = data$Age,
    Type = gsub("_DNA_weight", "", dna_col_name), 
    mean_biomass = rowMeans(sim_results, na.rm = TRUE),
    median_biomass = apply(sim_results, 1, median, na.rm = TRUE),
    q1_biomass = apply(sim_results, 1, quantile, probs = 0.25, na.rm = TRUE),
    q3_biomass = apply(sim_results, 1, quantile, probs = 0.75, na.rm = TRUE)
  )
}

# --- B. Sobol Sensitivity Analysis ---
run_sobol_analysis <- function(data, dna_col_name, organism_name, conversion_data, n_sobol) {
  factors <- conversion_data %>% filter(Organism == organism_name)
  if(nrow(factors) == 0 || !(dna_col_name %in% names(data))) return(NULL)
  
  dna_pool <- as.numeric(unlist(factors %>% select(starts_with("DNA_value")))); dna_pool <- dna_pool[!is.na(dna_pool)]
  bio_pool <- as.numeric(unlist(factors %>% select(starts_with("Biomass_value")))); bio_pool <- bio_pool[!is.na(bio_pool)]
  mean_dna_weight <- mean(as.numeric(as.character(data[[dna_col_name]])), na.rm=TRUE)
  if(is.na(mean_dna_weight) || mean_dna_weight == 0) return(NULL)
  
  d_min <- min(dna_pool); d_max <- max(dna_pool); b_min <- min(bio_pool); b_max <- max(bio_pool)
  
  struct_base <- 1.0
  if(grepl("Bacteria", organism_name)) struct_base <- PARAM_STRUCT_BACTERIA
  if(grepl("Algae", organism_name)) struct_base <- PARAM_STRUCT_ALGAE
  if(grepl("Woody", organism_name)) struct_base <- PARAM_STRUCT_WOODY
  if(grepl("Plant|Viridiplantae", dna_col_name)) struct_base <- PARAM_STRUCT_PLANT
  
  ext_min <- max(1.0, PARAM_EXTRACTION_EFFICIENCY - 1.0)
  ext_max <- PARAM_EXTRACTION_EFFICIENCY + 1.0
  
  X1 <- data.frame(P_DNA = runif(n_sobol, d_min, d_max), P_Bio = runif(n_sobol, b_min, b_max), P_Ext = runif(n_sobol, ext_min, ext_max), P_Str = runif(n_sobol, struct_base*0.8, struct_base*1.2))
  X2 <- data.frame(P_DNA = runif(n_sobol, d_min, d_max), P_Bio = runif(n_sobol, b_min, b_max), P_Ext = runif(n_sobol, ext_min, ext_max), P_Str = runif(n_sobol, struct_base*0.8, struct_base*1.2))
  
  model_fun <- function(X) { return( (mean_dna_weight / X[,1]) * X[,2] * X[,3] * X[,4] ) }
  x <- sensitivity::soboljansen(model = model_fun, X1 = X1, X2 = X2, nboot = 100)
  
  data.frame(Parameter = c("Genome Size", "Cell Biomass", "Extraction Eff.", "Structural Fact."), ST = x$T$original, Organism = organism_name)
}

# --- C. Perturbation Analysis (Corrected with different ranges) ---
run_perturbation_analysis <- function(data, dna_col_name, organism_name, conversion_data, n_sims) {
  res_base <- run_monte_carlo(data, dna_col_name, organism_name, conversion_data, n_sims)
  # ROBUSTNESS: Use Median to avoid outlier bias
  val_base <- median(res_base$median_biomass, na.rm=TRUE) 
  results <- list()
  
  # 1. Extraction (Methodological = Small Range)
  for(p in PERT_STEPS_METHOD) {
    res_p <- run_monte_carlo(data, dna_col_name, organism_name, conversion_data, n_sims, ext_factor = PARAM_EXTRACTION_EFFICIENCY * p)
    val_p <- median(res_p$median_biomass, na.rm=TRUE)
    results[[length(results)+1]] <- data.frame(Taxon=organism_name, Parameter="Extraction", Variation=paste0(p*100, "%"), Change_Pct = ((val_p - val_base)/val_base)*100)
  }
  
  # 2. Structure (Methodological = Small Range)
  base_structs <- list(bac=PARAM_STRUCT_BACTERIA, alg=PARAM_STRUCT_ALGAE, woody=PARAM_STRUCT_WOODY, plant=PARAM_STRUCT_PLANT)
  for(p in PERT_STEPS_METHOD) {
    mod_structs <- lapply(base_structs, function(x) x * p)
    res_p <- run_monte_carlo(data, dna_col_name, organism_name, conversion_data, n_sims, struct_factors = mod_structs)
    val_p <- median(res_p$median_biomass, na.rm=TRUE)
    results[[length(results)+1]] <- data.frame(Taxon=organism_name, Parameter="Structure", Variation=paste0(p*100, "%"), Change_Pct = ((val_p - val_base)/val_base)*100)
  }
  
  # 3. Genome Size (Biological = Large Range)
  for(p in PERT_STEPS_BIO) {
    res_p <- run_monte_carlo(data, dna_col_name, organism_name, conversion_data, n_sims, dna_scaling = p)
    val_p <- median(res_p$median_biomass, na.rm=TRUE)
    results[[length(results)+1]] <- data.frame(Taxon=organism_name, Parameter="Genome Size", Variation=paste0(p*100, "%"), Change_Pct = ((val_p - val_base)/val_base)*100)
  }
  
  # 4. Cell Biomass (Biological = Large Range)
  for(p in PERT_STEPS_BIO) {
    res_p <- run_monte_carlo(data, dna_col_name, organism_name, conversion_data, n_sims, bio_scaling = p)
    val_p <- median(res_p$median_biomass, na.rm=TRUE)
    results[[length(results)+1]] <- data.frame(Taxon=organism_name, Parameter="Cell Biomass", Variation=paste0(p*100, "%"), Change_Pct = ((val_p - val_base)/val_base)*100)
  }
  
  bind_rows(results)
}


# -----------------------------------------------------------------------------
# 3. MAIN EXECUTION LOOP
# -----------------------------------------------------------------------------

map_list <- list(
  list(col="Bacteria_aquatic_DNA_weight", org="Bacteria"),
  list(col="Bacteria_terrestrial_DNA_weight", org="Bacteria"),

  list(col="Archaea_aquatic_DNA_weight", org="Archaea"),
  list(col="Archaea_terrestrial_DNA_weight", org="Archaea"),

  list(col="Fungi_aquatic_DNA_weight", org="Fungi"),
  list(col="Fungi_terrestrial_DNA_weight", org="Fungi"),

  list(col="Metazoa_aquatic_DNA_weight", org="Metazoa"),
  list(col="Metazoa_terrestrial_DNA_weight", org="Metazoa"),

  list(col="Viridiplantae_woody_DNA_weight", org="Non_woody plant"),
  list(col="Viridiplantae_non_woody_DNA_weight", org="Non_woody plant"),
  list(col="Viridiplantae_aquatic_DNA_weight", org="Aquatic plant"),

  list(col="Viruses_DNA_weight", org="Viruses"),

  list(col="Aquatic_algae_DNA_weight", org="Algae")
)

all_cores_biomass <- list() 
plot_list_cores <- list() 
sobol_results_list <- list()
perturbation_results_list <- list()

for (core_id in core_list) {
  print(paste(">>> PROCESSING CORE", core_id))
  
  core_data <- DNA_global %>% filter(Number == core_id)
  names(core_data) <- names(core_data) %>% gsub("_Start_wt", "_DNA_weight", .) %>% gsub("_percentage", "", .)
  
  if(nrow(core_data) == 0) next  
  
  lake_name <- core_data$Lake[1]  
  lake_name_clean <- gsub("[^A-Za-z0-9]", "_", lake_name)  
  
  core_out_path <- file.path(PATH_OUTPUT_DIR, lake_name_clean)
  if(!dir.exists(core_out_path)) dir.create(core_out_path)
  
  mc_file <- file.path(core_out_path, paste0("MC_Results_", lake_name_clean, ".csv"))
  
  if(nrow(core_data) > 0) {
    mc_res <- list()
    for(m in map_list) {
      if(m$col %in% names(core_data)) {
        
        # 1. Main Simulation
        res <- run_monte_carlo(core_data, m$col, m$org, values_per_cell_global, N_SIMS_MC)
        if(!is.null(res)) mc_res[[length(mc_res)+1]] <- res
        
        # 2. Sobol Analysis
        sob <- run_sobol_analysis(core_data, m$col, m$org, values_per_cell_global, N_SOBOL)
        if(!is.null(sob)) {
          sob$Core <- core_id
          sobol_results_list[[length(sobol_results_list)+1]] <- sob
        }
        
        # 3. Perturbation Analysis
        # Only run for one representative core (e.g. Core 1) to save time
        if(core_id == core_list[1]) { 
          pert <- run_perturbation_analysis(core_data, m$col, m$org, values_per_cell_global, 2000) 
          if(!is.null(pert)) perturbation_results_list[[length(perturbation_results_list)+1]] <- pert
        }
      }
    }
    
    if(length(mc_res) > 0) {
      df_taxa <- bind_rows(mc_res)
      write.csv(df_taxa, mc_file, row.names=FALSE)
      df_taxa$Core <- core_id
      all_cores_biomass[[length(all_cores_biomass)+1]] <- df_taxa
      
      # Plotting Preparation
      df_final <- df_taxa %>%
        mutate(C_Factor = case_when(grepl("Bacteria", Type) ~ 0.52, grepl("Viruses", Type) ~ 0.68, grepl("Archaea", Type) ~ 0.50, grepl("woody", Type) ~ 0.48, TRUE ~ 0.50)) %>%
        mutate(Age_Group = round(Age, 4)) %>%
        group_by(Age_Group) %>%
        summarise(Age = mean(Age), Sim_Med = sum(median_biomass * C_Factor, na.rm=T)*100, Sim_Q1 = sum(q1_biomass * C_Factor, na.rm=T)*100, Sim_Q3 = sum(q3_biomass * C_Factor, na.rm=T)*100) %>% ungroup()
      
      df_meas <- DNA_global %>% filter(Number == core_id) %>% select(Age, TOC) %>% rename(Measured_TOC=TOC) %>% filter(!is.na(Measured_TOC)) %>%
        mutate(Meas_Min = Measured_TOC * 0.95, Meas_Max = Measured_TOC * 1.05)
      
      # Stats
      df_test <- inner_join(df_final, df_meas, by="Age")
      if(nrow(df_test) > 2) {
        corr_res <- cor.test(df_test$Sim_Med, df_test$Measured_TOC, method = "spearman")
        wilcox_res <- wilcox.test(df_test$Sim_Med, df_test$Measured_TOC, paired = TRUE)
        median_ratio <- median(df_test$Sim_Med / df_test$Measured_TOC, na.rm=T)
        mag_label <- ifelse(wilcox_res$p.value > 0.05, "Consistent", "Different")
        subtitle_text <- paste0("Trend: R=", round(corr_res$estimate, 2), " (p=", format.pval(corr_res$p.value, digits=2), ") | Mag: Ratio ", round(median_ratio, 2), "x (p=", format.pval(wilcox_res$p.value, digits=2), " ", mag_label, ")")
      } else { subtitle_text <- "Insufficient data" }
      
      # Plot
      p_corr <- ggplot() +
        geom_ribbon(data=df_meas, aes(x=Age, ymin=Meas_Min, ymax=Meas_Max, fill="Measured TOC (±5%)"), alpha=0.25) +
        geom_smooth(data = df_meas, aes(x = Age, y = Measured_TOC, color = "Measured TOC (±5%)"), method = "loess",span=0.5, se = FALSE) +
        geom_ribbon(data=df_final, aes(x=Age, ymin=Sim_Q1, ymax=Sim_Q3, fill="Simulated Biomass (IQR)"), alpha=0.15) +
        geom_smooth(data = df_final, aes(x = Age, y = Sim_Med, color = "Simulated Biomass (IQR)"), method = "loess",span=0.5, se = FALSE) +
        scale_color_manual(name="", values=c("Measured TOC (±5%)"="#5d2c04", "Simulated Biomass (IQR)"="#BC8F8F")) +
        scale_fill_manual(name="", values=c("Measured TOC (±5%)"="#5d2c04", "Simulated Biomass (IQR)"="#BC8F8F")) +
        annotation_logticks(sides="l") +
        theme_bw(base_size = 8) + theme(legend.position="top", panel.grid = element_blank()) +
        labs(title=paste("", lake_name), subtitle=subtitle_text, y="Weight percentage (wt%)", x="Age (k yrs)")
      
      ggsave(file.path(core_out_path, "Final_Comparison_grid.png"), p_corr, width=10, height=6)
      plot_list_cores[[as.character(core_id)]] <- p_corr
    
      df_habitat <- df_taxa %>%
        filter(Core == core_id) %>%  
        mutate(
          Habitat = case_when(
            grepl("aquatic|algae", Type, ignore.case = TRUE) ~ "Aquatic",
            grepl("terrestrial|woody|non_woody", Type, ignore.case = TRUE) ~ "Terrestrial",
            TRUE ~ NA_character_
          ),
          C_Factor = case_when(
            grepl("Bacteria", Type) ~ 0.52,
            grepl("Viruses", Type) ~ 0.68,
            grepl("Archaea", Type) ~ 0.50,
            grepl("woody", Type) ~ 0.48,
            TRUE ~ 0.50
          ),
          Age_Group = round(Age, 4)
        ) %>%
        filter(!is.na(Habitat))
      
      df_habitat_sum <- df_habitat %>%
        group_by(Age_Group, Habitat) %>%
        summarise(
          Med = sum(median_biomass * C_Factor, na.rm = TRUE),
          Q1 = sum(q1_biomass * C_Factor, na.rm = TRUE),
          Q3 = sum(q3_biomass * C_Factor, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        group_by(Age_Group) %>%
        mutate(
          Med_ral = (Med / sum(Med)) * 100,
          Q1_ral = (Q1 / sum(Q1)) * 100,
          Q3_ral = (Q3 / sum(Q3)) * 100
        ) %>%
        ungroup()
      
      
      p_habitat <- ggplot(df_habitat_sum,
                          aes(x = Age_Group,
                              y = Med_ral,
                              fill = Habitat)) +
        geom_area(position = "fill", alpha = 0.9) +
        scale_y_continuous(
          labels = scales::percent,
          name = "Percentage of OC DNA_projected weight percentage (%)"
        ) +
        scale_fill_manual(
          values = c(Aquatic = "#1f77b4",
                     Terrestrial = "#2ca02c")
        ) +
        theme_bw(base_size = 8) +
        theme(legend.position = "top", ) +
        labs(title = paste("", lake_name),
             x = "Age (kyrs)")
      
      ggsave(
        file.path(core_out_path, paste0("Core_", lake_name, "_Habitat.png")),
        p_habitat,
        width = 10, height = 6, dpi = 300
      )
      
      if (!exists("p_habitat_list")) p_habitat_list <- list()
      p_habitat_list[[as.character(core_id)]] <- p_habitat
      
      }
  }
}



# =============================================================================
# 4. MERGE MONTE CARLO BIOMASS RESULTS WITH ORIGINAL DATA
# =============================================================================

Processed_DNA <- read_csv(PATH_DATA_INPUT, show_col_types = FALSE)

mc_all <- bind_rows(all_cores_biomass)


mc_all <- mc_all %>%
  mutate(
    Number = Core 
  )

mc_wide <- mc_all %>%
  pivot_longer(
    cols = c(mean_biomass, median_biomass, q1_biomass, q3_biomass),
    names_to = "stat",
    values_to = "value"
  ) %>%
  unite("Type_stat", Type, stat, sep = "_") %>%
  pivot_wider(
    names_from = Type_stat,
    values_from = value
  )
Processed_DNA_with_MC <- Processed_DNA %>%
  left_join(mc_wide, by = c("Number", "Age"))


output_file <- file.path(dirname(PATH_DATA_INPUT), "Processed_DNA_with_MC.csv")

PC_new <- read_csv(
  "D:/document/DOC/AWI/Lake organic/data/Plat_PCARDA_with_scores_60.csv",
  show_col_types = FALSE
)
colnames(PC_new)
Processed_DNA_with_MC$PC1 <- PC_new$PC1
Processed_DNA_with_MC$PC2 <- PC_new$PC2

output_file <- file.path(
  dirname(PATH_DATA_INPUT),
  "Processed_DNA_with_MC_biomass.csv"
)

write_csv(Processed_DNA_with_MC, output_file)


# -----------------------------------------------------------------------------
# 4. FINAL PLOTS
# -----------------------------------------------------------------------------


# A. CORE GRID
if(length(plot_list_cores) > 0) {
  print(">>> CREATING GRID PLOT...")
  final_grid <- wrap_plots(plot_list_cores, ncol = 3) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Comparison of Simulated Biomass vs. Measured TOC",
      subtitle = "Trend: Pearson Correlation | Magnitude: Median Ratio & Wilcoxon Test",
      theme = theme(
        legend.position = "top",
        legend.direction = "horizontal",
        plot.title = element_text(size = 18, face = "bold")
      )
    )
  ggsave(file.path(PATH_OUTPUT_DIR, "All_Cores_Comparison_Grid.png"), final_grid, width = 18, height = 12, dpi=300)
}

#habitat
# B. HABITAT GRID
if(length(p_habitat_list) > 0) {  
  print(">>> CREATING HABITAT GRID PLOT...")
  
  final_habitat_grid <- wrap_plots(p_habitat_list, ncol = 3) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Relative Biomass per Habitat Across All Cores",
      subtitle = "Aquatic vs Terrestrial share over Age",
      theme = theme(
        legend.position = "bottom",
        plot.title = element_text(size = 18, face = "bold"))
    )
  
  ggsave(file.path(PATH_OUTPUT_DIR, "All_Cores_Habitat_Grid.png"), final_habitat_grid, width = 18, height = 12, dpi = 300)
}



# B. SOBOL PLOT

if(length(sobol_results_list) > 0) {
  print(">>> CREATING SOBOL PLOT...")
  df_sobol <- bind_rows(sobol_results_list)
  df_sobol <- df_sobol %>%
    mutate(Organism = ifelse(Organism == "Algae", "Algae/protists", Organism))
  p_sobol <- ggplot(df_sobol, aes(x = Parameter, y = ST, fill = Parameter)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ Organism, scales = "free_y") +
    theme_bw(base_size = 14) +
    labs(title = "Sobol Sensitivity Analysis", subtitle = "Which parameter drives uncertainty? (Variance Contribution)", y = "Total Sensitivity Index (ST)", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(PATH_OUTPUT_DIR, "Sobol_Sensitivity_Indices.png"), p_sobol, width = 14, height = 10, dpi=300)
}

# C. PERTURBATION PLOT
if(length(perturbation_results_list) > 0) {
  print(">>> CREATING PERTURBATION ANALYSIS PLOT...")
  df_pert <- bind_rows(perturbation_results_list)
  df_pert <- df_pert %>%
    mutate(Taxon = ifelse(Taxon == "Algae", "Algae/protists", Taxon))
  p_pert <- ggplot(df_pert, aes(x = Variation, y = Change_Pct, fill = Parameter)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~ Taxon, scales="free_y") +
    theme_bw(base_size = 14) +
    labs(title = "Perturbation Analysis (Stress Test)", subtitle = "Impact of changing assumptions on final Median Biomass", y = "Change in Median Biomass (%)", x = "Parameter Variation") +
    geom_hline(yintercept=0, color="black") +
    theme(legend.position="bottom")
  ggsave(file.path(PATH_OUTPUT_DIR, "Perturbation_Analysis_Bars.png"), p_pert, width = 16, height = 12, dpi=300)
  ggsave(file.path(PATH_OUTPUT_DIR, "Perturbation_Analysis_Bars.PDF"), p_pert, width = 16, height = 12, dpi=300)
}

#D. GLOBAL BARPLOT
if(length(all_cores_biomass) > 0) {
  print(">>> CREATING GLOBAL AGGREGATION PLOT...")
  df_biomass_global <- bind_rows(all_cores_biomass) %>% mutate(Taxon = gsub("_DNA_weight", "", Type), Taxon = gsub("_", " ", Taxon), Value = median_biomass * 100, Source = "Biomass weight percentage") %>% filter(Taxon %in% taxon_order) 
  df_dna_global <- DNA_global %>% select(ends_with("_DNA_weight")) %>% pivot_longer(cols = everything(), names_to = "Type", values_to = "Value") %>% mutate(Taxon = gsub("_DNA_weight", "", Type), Taxon = gsub("_", " ", Taxon), Value = Value * 100, Source = "DNA weight percentage") %>% filter(Taxon %in% taxon_order) 
  
  df_combined_global <- bind_rows(df_dna_global %>% select(Taxon, Value, Source), df_biomass_global %>% select(Taxon, Value, Source)) %>% group_by(Taxon, Source) %>% summarise(median = median(Value, na.rm=TRUE), q1 = quantile(Value, 0.25, na.rm=TRUE), q3 = quantile(Value, 0.75, na.rm=TRUE), .groups = "drop") %>% mutate(Taxon = factor(Taxon, levels = taxon_order))
  df_comp_global <- df_biomass_global %>% group_by(Core, Age) %>% mutate(Total_Bio = sum(Value, na.rm=TRUE)) %>% ungroup() %>% mutate(Percentage = (Value / Total_Bio) * 100) %>% group_by(Taxon) %>% summarise(median = median(Percentage, na.rm=TRUE), q1 = quantile(Percentage, 0.25, na.rm=TRUE), q3 = quantile(Percentage, 0.75, na.rm=TRUE), .groups = "drop") %>% mutate(Taxon = factor(Taxon, levels = taxon_order))
  
  pseudo_log_trans <- function(base = 10, sigma = 1) { trans <- function(x) sign(x) * log10(1 + abs(x) / sigma); inv <- function(x) sign(x) * sigma * (10^abs(x) - 1); trans_new("pseudo_log", trans, inv, domain = c(0, Inf)) }
  
  p_global_top <- ggplot(df_combined_global, aes(x = Taxon, y = median, fill = Taxon, color = Taxon)) + geom_col(aes(fill = ifelse(Source == "Biomass weight percentage", as.character(Taxon), NA), color = as.character(Taxon), group = Source, linewidth = ifelse(Source == "DNA weight percentage", 1.0, 0.3)), position = position_dodge(width = 0.7), width = 0.6, alpha=0.8) + geom_errorbar(aes(ymin = q1, ymax = q3, color = Taxon, group = Source), position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.8) + scale_fill_manual(values = setNames(color_values_alpha, taxon_order), na.value = "transparent") + scale_color_manual(values = setNames(color_values, taxon_order), guide = "none") + scale_linewidth_identity() + scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-6), breaks = c(1e-6, 1e-4, 1e-2, 1, 100), labels = trans_format("log10", math_format(10^.x))) + theme_bw(base_size = 16) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + labs(y = "Weight % (log10 scale)", title = "Global Aggregation")
  
  p_global_bot <- ggplot(df_comp_global, aes(x = Taxon, y = median, fill = Taxon, color = Taxon)) + geom_col(width = 0.6, linewidth = 0.8) + geom_errorbar(aes(ymin = q1, ymax = q3,color = Taxon), width = 0.2, linewidth = 0.8) + coord_cartesian(ylim = c(0, max(df_comp_global$q3, na.rm=TRUE)*1.1)) + scale_fill_manual(values = color_values_alpha) + scale_color_manual(values = color_values) + theme_bw(base_size = 16) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) + labs(y = "Biomass Share (%)", x = "")
  
  p_global_combined <- p_global_top / p_global_bot +
    plot_annotation(
      tag_levels = 'A',
      theme = theme(
        plot.tag = element_text(size = 20, face = "bold")
      )
    )
  
  ggsave(file.path(PATH_OUTPUT_DIR, "Global_Figure_4_Final_ExtractionCorrected.png"), p_global_combined, width = 14, height = 16, dpi=300)
  ggsave(file.path(PATH_OUTPUT_DIR, "Global_Figure_4_Final_ExtractionCorrected.pdf"), p_global_combined, width = 14, height = 16, dpi=300)
}

ratio_df <- df_combined_global %>%
  select(Taxon, Source, median) %>%
  pivot_wider(
    names_from = Source,
    values_from = median
  ) %>%
  mutate(
    Biomass_to_DNA_ratio = `Biomass weight percentage` /
      `DNA weight percentage`
  )

ratio_df
# =============================================================================
# 7. DIAGNOSTIC PLOT: PARAMETER SPACE
# =============================================================================
if(length(simulation_points_cache) > 0) {
  print(">>> CREATING PARAMETER SPACE PLOT...")
  df_diag_real <- if(is.data.frame(simulation_points_cache)) simulation_points_cache else bind_rows(simulation_points_cache)
  df_diag_real <- df_diag_real %>%
    mutate(Group = ifelse(Group == "Algae", "Algae/protists", Group))
  p_param_space <- ggplot(df_diag_real, aes(x = DNA_Factor, y = Cell_Biomass)) +
    geom_point(alpha = 0.2, size = 0.8, color = "#5c6d7e") +
    facet_wrap(~ Group, scales = "free") +
    geom_density_2d(data = subset(df_diag_real, Group == "Bacteria"), color = "#d62728", size = 0.6, bins = 6) +
    labs(title = "Parameter Space: 4-Way Logic Applied", subtitle = "Visualizing the Beta-Distributions used in MC Simulation", x = "DNA per Cell [g]", y = "Biomass C per Cell [g]") +
    theme_bw(base_size = 14) + theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold")) +
    scale_x_continuous(labels = scales::scientific) + scale_y_continuous(labels = scales::scientific)
  ggsave(file.path(PATH_OUTPUT_DIR, "Parameter_Space_4Way_LINEAR.png"), p_param_space, width = 16, height = 12, dpi=300)
}

print("ALL ANALYSES FINISHED SUCCESSFULLY!")
