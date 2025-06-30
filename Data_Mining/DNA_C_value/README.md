
# Selenium C-values Crawler & scraper

This Python script uses Selenium to scrape plant C-value data from [Kew’s C-values database](https://cvalues.science.kew.org/search).

## Features
- Opens the website using Safari WebDriver.
- Accepts cookies automatically.
- Selects the "Family" checkbox filter.
- Clicks the search button.
- Iterates through paginated search results (pages 55 to 123 by default).
- Extracts tabular data from each page.
- Saves the combined data as an Excel file (`C_values_Viridiplantae.xlsx`).
- Includes random delays and debug screenshots for stability and troubleshooting.

## Requirements
- Safari browser and Safari WebDriver installed and enabled.
- Python packages: `selenium`, `pandas`.

## Usage
1. Install dependencies if needed:
   ```bash
   pip install selenium pandas
   ```
2. Run the script:
   ```bash
   python cvalues_selenium_scraper.py
   ```
3. Data will be saved as `C_values_Viridiplantae.xlsx` in the current folder.
4. Debug screenshots are saved throughout to help monitor progress.
