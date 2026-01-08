import time
import re
import pandas as pd
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# --- CONFIGURATION ---
input_csv = "/Users/josefineweiss/Desktop/Zijuan Manuscript/Master_Bacteria_Raw_species.csv"  # Your input file
output_csv = "MicrobeAtlas_Full_Results_ALLBACTERIA.csv"         # Where the results will be saved
delimiter = ";"                                              # CSV delimiter (; or ,)

# --- BROWSER SETUP ---
options = webdriver.SafariOptions()
# Safari support for "--headless" varies; if it doesn't work, the window will remain visible.
options.add_argument("--headless") 

print("Starting Safari in background...")
driver = webdriver.Safari(options=options) 
driver.maximize_window()

# --- LOAD DATA ---
print(f"Loading file: {input_csv}")
try:
    # Attempt to load with the specified semicolon delimiter
    df = pd.read_csv(input_csv, delimiter=delimiter)
except Exception:
    print("Error loading with semicolon. Checking delimiter! Trying comma...")
    df = pd.read_csv(input_csv, delimiter=",")

print(f"--- START: Processing {len(df)} bacteria ---")

# --- MAIN LOOP ---
for index, row in df.iterrows():
    # Take the FIRST COLUMN (iloc[0]) as the species name
    species = str(row.iloc[0]).strip()
    
    # Skip empty rows
    if not species or species.lower() == "nan": 
        continue

    print(f"\n[{index+1}/{len(df)}] Processing: {species}")

    try:
        # 1. Open home page
        driver.get("https://microbeatlas.org/")
        
        # 2. Search
        search_box = WebDriverWait(driver, 8).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR, ".tagify__input"))
        )
        search_box.click()
        search_box.send_keys(species)
        time.sleep(1.2) # Wait for Tagify/Autocomplete to respond
        search_box.send_keys(Keys.RETURN)

        # 3. Wait for result list
        try:
            results = WebDriverWait(driver, 8).until(
                EC.presence_of_all_elements_located((By.CLASS_NAME, "data_row_record"))
            )
        except Exception:
            results = [] # Timeout -> no results found

        if len(results) > 0:
            # 4. Click the first match
            first_row = results[0]
            link = first_row.find_element(By.CSS_SELECTOR, "a.customTooltip")
            driver.execute_script("arguments[0].click();", link)
            
            # 5. Wait for statistics (Colored circles)
            WebDriverWait(driver, 15).until(
                EC.presence_of_element_located((By.CLASS_NAME, "sample_detailed_col"))
            )
            time.sleep(2.5) # Wait for animations to finish
            
            # 6. PARSE HTML & CREATE DYNAMIC COLUMNS
            soup = BeautifulSoup(driver.page_source, 'html.parser')
            # Regex to find: Category Name: Numbers (Percentage%)
            pattern = re.compile(r"(.*?):\s*[\d,.]+\s*(?:samples)?\s*\(([\d.]+)%\)", re.IGNORECASE)
            
            columns = soup.find_all("div", class_="sample_detailed_col")
            
            # Temporary storage for current species data
            current_species_data = {}

            for col in columns:
                # Main category (Header)
                header_div = col.find("div", class_="d-flex")
                if header_div:
                    header_text = header_div.get_text(strip=True)
                    match = pattern.search(header_text)
                    if match:
                        main_cat = match.group(1).strip().title() # e.g., "Aquatic"
                        main_val = float(match.group(2))
                        
                        # Column name: Total_Aquatic
                        current_species_data[f"Total_{main_cat}"] = main_val

                        # Subcategories (Lines)
                        sub_lines = col.find_all("div", class_="sub_env_line")
                        for sub in sub_lines:
                            sub_text = sub.get_text(strip=True)
                            sub_match = pattern.search(sub_text)
                            if sub_match:
                                sub_name = sub_match.group(1).strip().title() # e.g., "Marine"
                                sub_val = float(sub_match.group(2))
                                
                                # Column name: Aquatic_Marine (Replacing spaces with underscores)
                                clean_main = main_cat.replace(" ", "_")
                                clean_sub = sub_name.replace(" ", "_")
                                col_name = f"{clean_main}_{clean_sub}"
                                
                                current_species_data[col_name] = sub_val

            # 7. UPDATE DATAFRAME
            for col_name, val in current_species_data.items():
                # Initialize column if it doesn't exist
                if col_name not in df.columns:
                    df[col_name] = 0.0
                
                # Assign value
                df.at[index, col_name] = val
            
            print(f"   -> Success! ({len(current_species_data)} categories loaded)")
            df.at[index, "MA_Status"] = "Found"

        else:
            print("   -> No results found.")
            df.at[index, "MA_Status"] = "No Hits"

    except Exception as e:
        print(f"   -> Error: {e}")
        df.at[index, "MA_Status"] = "Error"

    # 8. SAVE (Progressively after each bacterium)
    # This overwrites the file with the current state including new dynamic columns
    df.to_csv(output_csv, index=False, sep=";")
    
    # Brief pause to reduce server load/browser stress
    time.sleep(1)

driver.quit()
print(f"\nFINISHED! Results saved in: {output_csv}")