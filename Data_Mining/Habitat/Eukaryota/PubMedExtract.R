library(rentrez)
library(dplyr)
library(stringr)

# List of species
species_list <- Eukaryota_species_List$Species

# Keywords for Habitat_Keywords detection
keywords <- c("marine", "freshwater", "river", "pond", "soil", "terrestrial", "lake")

# Function to extract keywords from abstracts
extract_keywords <- function(text, keywords) {
  found_keywords <- keywords[sapply(keywords, function(k) grepl(k, text, ignore.case = TRUE))]
  if (length(found_keywords) > 0) {
    return(paste(found_keywords, collapse = ", "))  
  } else {
    return(NA)  
  }
}

# Function to fetch abstracts with random timeouts
get_abstracts_with_Habitat_Keywords <- function(species, max_results = 5) {
  query <- paste0('"', species, '"')
  
  # Random delay before sending request (between 0.5 and 2 seconds)
  Sys.sleep(runif(1, 0.5, 2))
  
  search_results <- tryCatch(
    entrez_search(db = "pubmed", term = query, retmax = max_results),
    error = function(e) return(NULL)  
  )
  
  if (is.null(search_results) || length(search_results$ids) == 0) {
    return(NULL)
  }
  
  pmid_list <- search_results$ids
  
  abstracts <- lapply(pmid_list, function(pmid) {
    Sys.sleep(runif(1, 0.5, 2))  # Random delay before fetching each abstract
    abs_text <- tryCatch(
      entrez_fetch(db = "pubmed", id = pmid, rettype = "abstract", retmode = "text"),
      error = function(e) return(NA)
    )
    
    found_keywords <- extract_keywords(abs_text, keywords)
    
    return(data.frame(
      Species = species, 
      PMID = pmid, 
      Abstract = abs_text, 
      Habitat_Keywords_Keywords = found_keywords,
      stringsAsFactors = FALSE
    ))
  })
  
  return(do.call(rbind, abstracts))
}

# Batch processing (200 species per batch)
batch_size <- 200
num_batches <- ceiling(length(species_list) / batch_size)
batch_results <- list()

for (i in seq_len(num_batches)) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(species_list))
  species_batch <- species_list[start_idx:end_idx]
  
  message(paste("Processing batch", i, "of", num_batches, "(", start_idx, "-", end_idx, ")"))
  
  batch_data <- lapply(species_batch, function(species) {
    message(paste("Fetching:", species))
    abstracts_df <- get_abstracts_with_Habitat_Keywords(species, max_results = 1)
    
    if (is.null(abstracts_df)) {
      message(paste("No results for", species))
      return(NULL)
    } else {
      message(paste("Results for", species, ":", nrow(abstracts_df)))
    }
    
    return(abstracts_df[!is.na(abstracts_df$Habitat_Keywords), ])
  })
  
  batch_results[[i]] <- do.call(rbind, batch_data)  
}

# Combine all results
final_results <- do.call(rbind, batch_results)

# Function for habitat assignment of Habitat_Keywords
habitat_categorization <- function(data) {
  data %>%
    mutate(
      category = case_when(
        # marine, freshwater & terrestrial existing → marine/freshwater/terrestrial
        str_detect(Habitat_Keywords, regex("marine", ignore_case = TRUE)) & 
          str_detect(Habitat_Keywords, regex("freshwater|pond|river|lake", ignore_case = TRUE)) & 
          str_detect(Habitat_Keywords, regex("terrestrial|soil", ignore_case = TRUE)) ~ "marine/freshwater/terrestrial",
        
        # freshwater & terrestrial existing → freshwater/terrestrial
        str_detect(Habitat_Keywords, regex("freshwater|pond|river|lake", ignore_case = TRUE)) & 
          str_detect(Habitat_Keywords, regex("terrestrial|soil", ignore_case = TRUE)) ~ "freshwater/terrestrial",
        
        # marine & freshwater existing → marine/freshwater
        str_detect(Habitat_Keywords, regex("marine", ignore_case = TRUE)) & 
          str_detect(Habitat_Keywords, regex("freshwater|pond|river|lake", ignore_case = TRUE)) ~ "marine/freshwater",
        
        # marine & terrestrial existing → marine/terrestrial
        str_detect(Habitat_Keywords, regex("marine", ignore_case = TRUE)) & 
          str_detect(Habitat_Keywords, regex("terrestrial|soil", ignore_case = TRUE)) ~ "marine/terrestrial",
        
        # only freshwater → freshwater
        str_detect(Habitat_Keywords, regex("freshwater|pond|river|lake", ignore_case = TRUE)) ~ "freshwater",
        
        # only terrestrial → terrestrial
        str_detect(Habitat_Keywords, regex("terrestrial|soil", ignore_case = TRUE)) ~ "terrestrial",
        
        # only marine → marine
        str_detect(Habitat_Keywords, regex("marine", ignore_case = TRUE)) ~ "marine",
        
        # If nothing fits → NA
        TRUE ~ NA_character_
      )
    )
}


# use function for dataset 
final_results_habitat <- habitat_categorization(final_results)

# show results
unique(final_results_habitat$category)
# save results
write.csv(final_results_habitat, "final_results_habitat.csv", row.names = FALSE)

message("Done! Results saved in 'species_Habitat_Keywords_results.csv'")

