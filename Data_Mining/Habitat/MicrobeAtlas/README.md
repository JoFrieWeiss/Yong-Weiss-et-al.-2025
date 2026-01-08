# MicrobeAtlas Environmental Data Scraper

This document contains the project overview, setup instructions, and the complete Python source code for the MicrobeAtlas scraper.

---

## 1. Description
This script automates the extraction of environmental distribution data from [MicrobeAtlas.org](https://microbeatlas.org/). It maps bacterial species names to their natural habitats (e.g., Aquatic, Terrestrial, Host-associated) and records their distribution percentages.



## 2. Features
* **Automated Search & Navigation:** Uses Selenium (Safari) to interact with the web interface.
* **Dynamic Column Creation:** Automatically detects and adds columns like `Total_Aquatic` or `Aquatic_Marine` based on the search results.
* **Progressive Auto-Save:** Writes to the CSV file after every processed species to ensure data safety.
* **Status Tracking:** Updates an `MA_Status` column with "Found", "No Hits", or "Error".

---

## 3. Prerequisites

### macOS Safari Setup
1.  Open **Safari**.
2.  Go to **Settings** (or Preferences) > **Advanced**.
3.  Check **"Show Develop menu in menu bar"**.
4.  In the **Develop** menu, click **"Allow Remote Automation"**.

### Python Libraries
Run the following command in your terminal:
```bash
pip install pandas beautifulsoup4 selenium
