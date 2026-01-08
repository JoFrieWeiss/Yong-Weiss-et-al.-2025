
# BacDive Database Scraper

This Python script scrapes bacterial strain data from the BacDive database and enriches an input CSV file with additional metadata.

## Features

- Reads an input CSV file containing bacterial species names.
- Searches BacDive for each species.
- Handles multiple search hits by selecting the first strain result.
- Extracts metadata such as:
  - BacDive URL for the strain
  - Isolation source categories
  - Optimal pH value
  - Sample counts from different environments (aquatic, soil, animal, plant)
- Updates the input CSV file with the new data.
- Prints progress and handles timeouts to avoid overloading the server.

## Requirements

- Python 3
- Packages: `re`, `time`, `pandas`, `requests`, `beautifulsoup4`

Install missing packages with:

```bash
pip install pandas requests beautifulsoup4
```

## Usage

1. Place your input CSV file in the same directory as the script. The CSV must have species names in the first column and use semicolon (`;`) as delimiter.

2. Set the filename in the script variable `csv_filename`.

3. Run the script:

```bash
python bacdive_scraping_v2.1.py
```

4. The script will create a new CSV file prefixed with `updated_` containing the additional scraped information.

## Notes

- The script pauses every 1000 requests (configurable by `batch`) for 5 seconds (`timeout`) to avoid hitting server rate limits.
- If no hits are found for a species, it logs this to the console.
- The script updates the CSV incrementally to save progress in case of interruption.

---

Feel free to modify batch size, timeout, or columns to extract based on your needs.


**Author:** Josefine Friederike Weiß
**Date:** 2025

