import zipfile
import os
import glob

def get_distribution(zip_path):
    """
    Extracts the Sequence Length Distribution module from a FastQC zip file.
    """
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            # Search for the fastqc_data.txt file within the zip archive
            data_filename = next((n for n in z.namelist() if n.endswith('fastqc_data.txt')), None)
            if not data_filename:
                print(f"  [!] fastqc_data.txt not found in: {zip_path}")
                return None

            with z.open(data_filename) as f:
                lines = f.read().decode('utf-8').splitlines()

        results = []
        is_in_block = False
        
        for line in lines:
            # Look for the start of the relevant module
            if ">>Sequence Length Distribution" in line:
                is_in_block = True
                continue
            
            # Look for the end of the module
            if is_in_block and line.startswith(">>END_MODULE"):
                break
            
            # Collect data (exclude headers starting with '#' and empty lines)
            if is_in_block and not line.startswith("#") and line.strip():
                results.append(line.strip())
        
        return results

    except Exception as e:
        print(f"  [!] Error processing {zip_path}: {e}")
        return None

# Define output filename
output_file = "length_distribution_results.txt"

# Search for all zip files in the current directory
zips = sorted(glob.glob("*.zip"))
print(f"Zip files found in the current directory: {len(zips)}")

if len(zips) == 0:
    print("ERROR: No .zip files found!")
else:
    with open(output_file, "w") as out:
        for zip_file in zips:
            print(f"Processing: {zip_file}...")
            data = get_distribution(zip_file)
            if data:
                out.write(f"Sample: {zip_file}\n")
                out.write("\n".join(data) + "\n")
                out.write("-" * 40 + "\n")
            else:
                print(f"  [?] No data extracted for {zip_file}")

    print(f"\nDone! Results written to: {output_file}")