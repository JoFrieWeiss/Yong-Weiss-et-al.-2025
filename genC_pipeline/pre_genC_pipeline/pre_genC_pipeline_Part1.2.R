### ILIRNEY DATA PROCESSING ###

# Load and filter data
Ilirney <- MergedData
Ilirney_filtered <- Ilirney %>% filter(ka > 0)

# Total reads per time point
Ilirney_summary <- Ilirney_filtered %>%
  group_by(ka) %>%
  summarise(total_reads = sum(taxonReads, na.rm = TRUE))

# Reads without assigned superkingdom
Ilirney_super_na <- Ilirney %>% filter(is.na(superkingdom) & ka > 0)
Ilirney_super_na_summary <- Ilirney_super_na %>%
  group_by(ka) %>%
  summarise(super_na_reads = sum(taxonReads, na.rm = TRUE))

# Subtract unassigned reads to get final totals
Ilirney_cleaned <- left_join(Ilirney_summary, Ilirney_super_na_summary, by = "ka") %>%
  mutate(super_na_reads = replace_na(super_na_reads, 0),
         reads_without_na = total_reads - super_na_reads)

# Assign habitats to species
habitat_reference <- read_excel("D:/document/DOC/AWI/Lake organic/data/final_results_eukaryota_excel.xlsx") %>%
  rename(species = Species)

eukaryota_Ilirney <- Ilirney %>%
  filter(superkingdom == "Eukaryota", rank == "species", ka > 0)

Ilirney_annotated <- left_join(eukaryota_Ilirney, habitat_reference, by = "species") %>%
  mutate(new_habitat = case_when(
    phylum == "Streptophyta" & !family %in% c("Butomaceae", "Characeae", "Closteriaceae", "Desmidiaceae",
                                              "Elatinaceae", "Equisetaceae", "Hydrocharitaceae", "Juncaceae",
                                              "Nymphaeaceae", "Lentibulariaceae", "Potamogetonaceae",
                                              "Typhaceae", "Juncaginaceae") ~ "terrestrial",
    order == "Lepidoptera" ~ "terrestrial",
    family == "Potamogetonaceae" ~ "freshwater",
    species == "Allium schoenoprasum" ~ "terrestrial",
    species %in% c("Equisetum hyemale", "Luzula luzuloides", "Sparganium erectum", "Triglochin maritima") ~ "freshwater/terrestrial",
    species == "Penaeus chinensis" ~ "marine",
    order == "Octopoda" ~ "marine",
    family == "Cordylidae" ~ "terrestrial",
    genus == "Cherax" ~ "freshwater",
    species %in% c("Boleophthalmus pectinirostris", "Salvelinus fontinalis") ~ "freshwater",
    TRUE ~ Habitat
  ))

### FUNGI PROCESSING ###

# All fungal reads
fungi_all <- Ilirney_filtered %>% filter(kingdom == "Fungi")
fungi_summary <- fungi_all %>%
  group_by(ka) %>%
  summarise(fungi_total = sum(taxonReads, na.rm = TRUE))

# Habitat-specific fungal reads
fungi_annotated <- Ilirney_annotated %>% filter(kingdom == "Fungi")

process_habitat_counts <- function(data, group_label) {
  grouped <- data %>%
    group_by(new_habitat, ka) %>%
    summarise(sumcount = sum(taxonReads), .groups = "drop") %>%
    pivot_wider(names_from = new_habitat, values_from = sumcount, values_fill = 0) %>%
    pivot_longer(-ka, names_to = "new_habitat", values_to = !!group_label)
  
  return(grouped)
}

fungi_by_habitat <- process_habitat_counts(fungi_annotated, "fungi_sum")

# Total fungi reads across habitats
fungi_habitat_total <- fungi_by_habitat %>%
  group_by(ka) %>%
  summarise(fungi_habitat_sum = sum(fungi_sum, na.rm = TRUE))

# Classify habitat groups
fungi_grouped <- left_join(Ilirney_summary, fungi_by_habitat, by = "ka") %>%
  mutate(habitat_group = case_when(
    new_habitat %in% c("freshwater", "marine", "marine/freshwater") ~ "aquatic",
    new_habitat == "terrestrial" ~ "terrestrial",
    new_habitat == "freshwater/terrestrial" ~ "freshwater_terrestrial",
    TRUE ~ "varied"
  )) %>%
  group_by(ka, habitat_group) %>%
  summarise(total_fungi = sum(fungi_sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = habitat_group, values_from = total_fungi, values_fill = 0)

# Merge all Fungi data
fungi_combined <- reduce(list(Ilirney_cleaned, fungi_summary, fungi_habitat_total, fungi_grouped),
                         full_join, by = "ka") %>%
  mutate(
    ter_aqu_total = rowSums(across(c(terrestrial, aquatic), ~ replace_na(.x, 0))),
    fungi_na_estimate = (fungi_total / reads_without_na) * super_na_reads,
    fungi_unclassified = (fungi_total - ter_aqu_total) + fungi_na_estimate,
    perc_fungi_aquatic = aquatic / total_reads,
    perc_fungi_terrestrial = terrestrial / total_reads,
    perc_fungi_unclassified = fungi_unclassified / total_reads,
    perc_fungi_total = fungi_total / total_reads
  )

# Save Fungi summary
write.csv(fungi_combined, "D:/document/DOC/AWI/Lake organic/data/Ilirney/Ilirney_Fungi_ter_aqu_fre_ter_percentage.csv", row.names = FALSE)

### VIRIDIPLANTAE PROCESSING ###

# Total Viridiplantae reads
viridi_all <- Ilirney_filtered %>% filter(kingdom == "Viridiplantae")
viridi_summary <- viridi_all %>%
  group_by(ka) %>%
  summarise(viridiplantae_total = sum(taxonReads, na.rm = TRUE))

# Habitat-specific Viridiplantae reads
viridi_annotated <- Ilirney_annotated %>% filter(kingdom == "Viridiplantae")

viridi_by_habitat <- process_habitat_counts(viridi_annotated, "viridi_sum")

# Total Viridiplantae reads across habitats
viridi_habitat_total <- viridi_by_habitat %>%
  group_by(ka) %>%
  summarise(viridi_habitat_sum = sum(viridi_sum, na.rm = TRUE))

# You can proceed with similar calculations for Viridiplantae as done for Fungi.
# Example:
# - classify habitat groups
# - calculate percentages
# - merge and save
