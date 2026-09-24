# 1. Load libraries
library(tidyverse)

# 2. Define main directory
base_path <- "/Users/josefineweiss/Desktop/Zijuan Manuscript/Read_length"
target_dirs <- c("Btoko", "Illirney", "Lama", "Lele", "Salmon", "Ulu")
search_paths <- file.path(base_path, target_dirs)

# 3. Find all TXT files
files <- list.files(path = search_paths, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
message("Files found: ", length(files))

# 4. Function to read the custom structure including folder name
read_fastqc_txt <- function(file_path) {
  # Extract the folder name (e.g., Btoko or Lama) from the path
  # We take the directory directly above the file
  folder_name <- basename(dirname(file_path))
  
  lines <- readLines(file_path)
  
  full_data <- list()
  current_sample <- NULL
  temp_rows <- list()
  
  for (line in lines) {
    if (startsWith(line, "Sample:")) {
      # Save previous sample, if present
      if (!is.null(current_sample) && length(temp_rows) > 0) {
        full_data[[current_sample]] <- bind_rows(temp_rows)
      }
      # Set new sample
      current_sample <- str_remove(line, "Sample: ") %>% str_trim()
      temp_rows <- list()
    } else if (grepl("^[0-9]", line)) {
      # Data row (Length \t Count)
      parts <- str_split(line, "\t")[[1]]
      if (length(parts) >= 2) {
        temp_rows[[length(temp_rows) + 1]] <- tibble(
          Group = folder_name,  # Store the folder name here
          Sample = current_sample,
          Length = parts[1],
          Count = as.numeric(parts[2])
        )
      }
    }
  }
  # Add the last sample of the file
  if (!is.null(current_sample) && length(temp_rows) > 0) {
    full_data[[current_sample]] <- bind_rows(temp_rows)
  }
  
  return(bind_rows(full_data))
}

# 5. Read and combine all files
raw_data <- files %>%
  map_dfr(read_fastqc_txt)

# 6. Calculate average (now grouped by Group and Sample)
results <- raw_data %>%
  separate(Length, into = c("Min_Len", "Max_Len"), sep = "-", fill = "right", convert = TRUE) %>%
  mutate(
    Effective_Length = ifelse(is.na(Max_Len), Min_Len, (Min_Len + Max_Len) / 2)
  ) %>%
  group_by(Group, Sample) %>% # Keep both
  summarise(
    Avg_Length = sum(Effective_Length * Count) / sum(Count),
    Total_Reads = sum(Count),
    .groups = "drop"
  ) %>%
  arrange(Group, Sample)

# 7. Display & save results
print(results)

output_path <- file.path(base_path, "average_sequence_lengths_with_Groups.csv")
write_csv2(results, output_path)

message("\nDone! The 'Group' column now contains the name of the original folder.")
message("File saved under: ", output_path)


# 1. Load libraries
library(tidyverse)

# 2. Define paths
base_path_json <- "/Users/josefineweiss/Desktop/Zijuan Manuscript/Reads_all"
target_dirs_json <- c("Btoko", "Illirney", "Lama", "Lele", "Salmon")
search_paths_json <- file.path(base_path_json, target_dirs_json)

# 3. Find all CSV files in these folders
csv_files <- list.files(path = search_paths_json, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
message("CSV files found in Reads_all: ", length(csv_files))

# 4. Read files and add folder name as group
json_csv_data <- csv_files %>%
  map_dfr(function(file_path) {
    # Extract folder name
    folder_name <- basename(dirname(file_path))
    
    # Read file
    read_csv(file_path, show_col_types = FALSE) %>%
      mutate(Group = folder_name) %>% # Add Group column
      select(Group, everything())     # Set Group as the first column
  })

# 5. Display result
print(json_csv_data)

# 6. Save
output_final <- file.path(base_path_json, "summary_json_parser_ALL_GROUPS.csv")
write_csv2(json_csv_data, output_final)

message("\nDone! The merged table is located at: ", output_final)