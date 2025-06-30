from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
import time
import random

# Webdriver setup
options = webdriver.SafariOptions()
driver = webdriver.Safari(options=options)

try:
    # Load webpage
    url = "https://cvalues.science.kew.org/search"
    driver.get(url)
    time.sleep(random.uniform(2, 5))  # Random wait time (2-5 seconds)

    # Save debug screenshot
    driver.save_screenshot("debug_step1_landing_page.png")

    # Accept cookies
    try:
        cookie_button = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//button[contains(text(), 'I Accept')]"))
        )
        cookie_button.click()
        print("Cookies accepted.")

        # Wait until cookie window fully disappears
        time.sleep(random.uniform(2, 4))  # Random wait time
        driver.save_screenshot("debug_step2_after_cookies.png")
    except Exception as e:
        print("No 'Accept Cookies' button found or error clicking it:", e)

    # Activate "Family" checkbox
    try:
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.XPATH, "//*[@id='familyCheckbox']"))
        )
        family_checkbox = driver.find_element(By.XPATH, "//*[@id='familyCheckbox']")
        driver.execute_script("arguments[0].scrollIntoView(true);", family_checkbox)  # Scroll to checkbox
        time.sleep(random.uniform(1, 3))  # Random pause for stability
        family_checkbox.click()
        print("Checkbox 'Family' clicked.")
        driver.save_screenshot("debug_step3_after_checkbox.png")
    except Exception as e:
        print("Checkbox 'Family' not found or error clicking it:", e)
        driver.save_screenshot("debug_no_family_checkbox.png")

    # Click on "Search" button
    try:
        search_button = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//input[@type='submit' and @value='Search']"))
        )
        search_button.click()
        print("Search button clicked.")
        driver.save_screenshot("debug_step4_after_search.png")
    except Exception as e:
        print("Search button not found or error clicking it:", e)
        driver.save_screenshot("debug_no_search_button.png")

    # Loop through all pages and extract data
    data = []
    page_number = 55  # Start at page 55
    max_pages = 123  # Maximum of 123 pages to scrape

    while page_number <= max_pages:
        # Wait for the table to appear
        WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.CSS_SELECTOR, "table.table.table-striped.search-results.table-sm"))
        )

        # Extract the table
        table = driver.find_element(By.CSS_SELECTOR, "table.table.table-striped.search-results.table-sm")
        rows = table.find_elements(By.TAG_NAME, "tr")

        # Extract data from rows
        for row in rows[1:]:  # Skip first row (header)
            cols = row.find_elements(By.TAG_NAME, "td")
            data.append([col.text for col in cols])

        # Save debug screenshot after data extraction
        driver.save_screenshot(f"debug_step5_page_{page_number}.png")

        # Random pause after data extraction
        time.sleep(random.uniform(3, 6))  # Random pause between 3 and 6 seconds

        # Click on the next page (based on data-page attribute)
        try:
            next_page_button = WebDriverWait(driver, 15).until(
                EC.element_to_be_clickable((By.XPATH, f"//a[@class='page-link cvalues-pagination' and @data-page='{page_number + 1}']"))
            )
            next_page_button.click()
            print(f"Moving to page {page_number + 1}...")
            page_number += 1

            # Random pause after page change
            time.sleep(random.uniform(5, 8))  # Random pause between 5 and 8 seconds

        except Exception as e:
            print("No more pages or error navigating:", e)
            break

    # Column names
    header = ["Family", "Genus", "Species", "Subspecies", "1C DNA Amount (pg)", "Original Reference"]

    # Load into DataFrame
    df = pd.DataFrame(data, columns=header)

    # Save to Excel
    df.to_excel("C_values_Viridiplantae.xlsx", index=False)
    print("Data successfully saved to 'C_values_Viridiplantae.xlsx'.")
finally:
    driver.quit()
