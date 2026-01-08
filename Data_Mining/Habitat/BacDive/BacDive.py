import re
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup

# Load csv to dataframe.
csv_filename = "YOURSPECIESLIST.csv"
df = pd.read_csv(csv_filename, delimiter=';')

# Define base urls - search and response.
bacdive_baseurl = "https://bacdive.dsmz.de/search?search="
strain_baseurl = "https://bacdive.dsmz.de/strain/"

# Set batch and timeout.
batch = 1000
timeout = 5

# Add new columns for data extraction.
df["bacdive_url"] = ""
df["isolation_categories"] = ""
df["pH_optimum"] = ""
df["aquatic_samples"] = ""
df["soil_samples"] = ""
df["animal_samples"] = ""
df["plant_samples"] = ""

# Initialize counters.
positive_count, negative_count = 0, 0

# Iterate over the dataframe.
for index, species in df.iloc[:, 0].items():
    # Concatenate search url from base and species.
    species_encoded = species.replace(" ", "+")
    url = bacdive_baseurl + species_encoded
    multiple = False

    # Fetch html content of the url.
    response = requests.get(url)
    html_content = response.text

    # Parse the HTML content.
    soup = BeautifulSoup(html_content, 'html.parser')

    # MULTIPLE HITS CASE - Replace URL with first hit.
    if "hits: " in html_content and not "hits: 0" in html_content:
        # Find list elements, take second link and create new url.
        list_results = soup.find_all("li", class_="searchresultrow1")
        second_list_element = list_results[1]
        first_link = second_list_element.find("a")
        strain_href = first_link.get("href")
        strain_href = strain_href.replace("/strain/", "")
        strain_first_hit_url = strain_baseurl + strain_href
        multiple = True
        # Need to create new soup with updated species url.
        url = strain_first_hit_url
        response = requests.get(url)
        html_content = response.text
        soup = BeautifulSoup(html_content, 'html.parser')

    # Open csv file and to append each record.
    output_filename = "updated_" + csv_filename
    with open(output_filename, 'a', newline='') as csvfile:

        # Check if response contains strain url or not.
        if strain_baseurl in response.url:
            if not "hits: " in html_content:
                # Print if multiple hits or just one.
                if multiple:
                    print(f"HIT (multiple): {url} - Multiple species found, taking first result.")
                else:
                    print(f"HIT (single): {url} - One species found.")
                # Update the dataframe with response url if yes.
                df.at[index, 'bacdive_url'] = response.url
                positive_count += 1
                all_categories = []

                # Debugging: Print HTML content
                #print("HTML content length:", len(html_content))

                # If table with categories exists, extract them and write to dataframe.
                isolation_sources_categories = soup.find('table', class_='detail-isol-categories')
                if isolation_sources_categories:
                    #print("Table found")
                    for td in isolation_sources_categories.find_all('td'):
                        if "#" in td.text:
                            all_categories.append(td.text.strip())
                    df.at[index, "isolation_categories"] = ", ".join(all_categories)

                # pH Optimum
                ph_table = soup.find('table', id="ph_table")
                if ph_table:
                    #print("Table found")
                    for row in ph_table.find_all('tr')[1:]:  # Skip header row
                        cells = row.find_all('td')
                        if len(cells) >= 4:  # Ensure enough columns
                            ph_type = cells[3].get_text(strip=True)
                            ph_value = cells[4].get_text(strip=True)
                            if ph_type == 'optimum':  # Only process if it's the optimum pH
                                # Use regex to extract the number directly into the DataFrame
                                df.at[index, "pH_optimum"] = re.search(r'(\d+(?:\.\d+)?)', ph_value).group(1)
                                break  # No need to check other rows once optimum is found

                # Extract Sample Count (By finding 'td' with class 'border')
                aquatic_samples_table = soup.find('table', class_='id_5')
                sample_types = ["Aquatic Samples", "Soil Samples", "Animal Samples", "Plant Samples"]

                if aquatic_samples_table:
                    #print("Aquatic samples table found")
                    for sample_type in sample_types:
                        for row in aquatic_samples_table.find_all('tr'):
                            for cell in row.find_all('td'):
                                if sample_type in cell.get_text(strip=True):
                                    count_cell = cell.find_next_sibling('td', class_='border')
                                    if count_cell:
                                        count = count_cell.get_text(strip=True)
                                        df.at[index, sample_type.lower().replace(" ", "_")] = count  # Create column name dynamically
                                    else:
                                        df.at[index, sample_type.lower().replace(" ", "_")] = None
                                        break  # Exit the inner loop after finding the sample

            # Save updated dataframe to CSV after each iteration.
            df.to_csv(csvfile, index=False, header=csvfile.tell()==0)

        # NO HITS CASE.
        else:
            # If really nothing is found.
            if "hits: 0" in html_content:
                print(f"no hit: {url} not found")
                negative_count += 1

        # Pause for some seconds every x requests.
        total_count = negative_count + positive_count
        if total_count % batch == 0:
            print(f"Pausing for {timeout} seconds after {total_count} responses.")
            time.sleep(timeout)

# Print final counts.
print("Found:" + str(positive_count), "Not found: " + str(negative_count))

# Save updated dataframe as csv.
df.to_csv("updated_" + csv_filename, index=False)

