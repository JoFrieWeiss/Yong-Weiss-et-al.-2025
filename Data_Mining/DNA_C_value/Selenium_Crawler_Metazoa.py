from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
import time
import random
import os

options = webdriver.SafariOptions()
driver = webdriver.Safari(options=options)

try:
    base_url = "https://www.genomesize.com/search.php"
    driver.get(base_url)
    time.sleep(random.uniform(2, 4))

    # All group links (text + URL)
    group_elements = driver.find_elements(By.XPATH, '//div[@id="search_by_type"]//a')
    group_links = [
        (el.text.strip().split(" (")[0].replace(":", "-").replace(" ", "_"), el.get_attribute("href"))
        for el in group_elements
    ]

    for group_name, href in group_links:
        print(f"===> Scraping group: {group_name}")
        driver.get(href)
        time.sleep(random.uniform(2, 4))

        data = []
        page_number = 1

        while True:
            print(f"{group_name} – page {page_number}")
            try:
                table = WebDriverWait(driver, 10).until(
                    EC.presence_of_element_located((By.ID, "query_results"))
                )
                rows = table.find_elements(By.TAG_NAME, "tr")
                for row in rows[1:]:
                    cols = row.find_elements(By.TAG_NAME, "td")
                    data.append([col.text.strip() for col in cols])

                # Check if there is another page
                next_page = page_number + 1
                next_button_xpath = f"//a[text()='{next_page}']"
                next_button = driver.find_elements(By.XPATH, next_button_xpath)

                if next_button:
                    driver.execute_script("arguments[0].scrollIntoView(true);", next_button[0])
                    time.sleep(random.uniform(1, 2))
                    next_button[0].click()
                    page_number += 1
                    time.sleep(random.uniform(3, 5))
                else:
                    print(f"{group_name} finished.")
                    break

            except Exception as e:
                print(f"Error in {group_name}, page {page_number}: {e}")
                break

        if data:
            df = pd.DataFrame(data)
            outdir = "GenomeSize_Groups"
            os.makedirs(outdir, exist_ok=True)
            df.to_excel(f"{outdir}/{group_name}.xlsx", index=False)
            print(f"→ Saved {group_name}.xlsx with {len(df)} rows.\n")

finally:
    driver.quit()
