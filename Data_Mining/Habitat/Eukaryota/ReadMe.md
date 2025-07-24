# Habitat Keywords Extraction Script

This R script searches PubMed abstracts for habitat-related keywords (e.g., marine, freshwater, terrestrial) for a list of species and assigns each species to habitat categories based on the found keywords.

## Features

- Fetches PubMed abstracts for each species (with random short delays to avoid server overload)  
- Searches abstracts for defined habitat keywords  
- Assigns habitat categories (marine, freshwater, terrestrial, or combinations) to species  
- Processes species in batches (200 species per batch)  
- Saves results to a CSV file  

## Requirements

R packages required:

- `rentrez` (for PubMed queries)  
- `dplyr` (for data manipulation)  
- `stringr` (for regex-based text search)  

## Usage

1. Define your species list in the variable `species_list` (e.g., from your dataframe).  
2. Run the script.  
3. Output: CSV file `final_results_habitat.csv` containing habitat keyword annotations and habitat category assignments.  

## Notes

The script uses random sleep intervals between requests to avoid overloading PubMed servers and to prevent IP blocking.

---

**Author:** Josefine Friederike Weiß
**Date:** 2025
