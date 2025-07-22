# HPC-Optimized Pipeline for Biomass Modeling
This repository contains a series of R scripts designed to process raw eDNA analysis data from sediment cores and convert it into biomass estimates. The pipeline is optimized for execution on a High-Performance Computing (HPC) cluster, leveraging parallel processing to significantly reduce computation time.

# Pipeline Overview
The workflow is divided into three main stages to achieve maximum efficiency on a cluster:

#Stage 1: Parallel Pre-processing (01_run_preprocessing.R)
This script processes the raw data (MergedData.RData) for a single sediment core.
It calculates the percentage abundance of various taxonomic groups (bacteria, fungi, plants, etc.) based on habitat and classification rules.
Goal: This computationally intensive step is run as a job array on the HPC cluster, allowing all cores to be processed simultaneously and independently.

#Stage 2: Combining Results (02_combine_preprocessing.R)
This script is executed once after all jobs from Stage 1 have completed successfully.
It gathers the intermediate results from all cores and merges them into a single, comprehensive CSV file.
It then calculates the DNA weights (..._Start_wt) based on concentration data from the Complex DNA information.csv file.
Goal: To create a clean, complete input file for the final biomass modeling stage.

#Stage 3: Parallel Biomass Modeling (03_run_biomass_model.R)
This script reads the combined file from Stage 2 and calculates the biomass permutations for a single sediment core.
It uses a cross_join to model all possible biomass values based on the uncertainties in DNA and biomass-per-cell values.
Goal: This intensive step is also run as a job array in parallel for all cores to generate the final results quickly.
Requirements
Ensure the following R packages are installed in your HPC environment:

# Run this in an R session to install any missing packages
install.packages(c("data.table", "readxl", "dplyr", "tidyr", "readr", "purrr", "stringr"))

Setup
Clone the Repository:
git clone [URL_OF_YOUR_REPOSITORY]
cd [YOUR_REPOSITORY_NAME]

Save Files: Ensure the three R scripts (01_..., 02_..., 03_...) and the bash script (run_hpc.sh) are in the main directory of your project.
Adjust Paths (IMPORTANT): Open each of the three R scripts and modify the BASE_DIR variable to match the directory structure of your HPC cluster. All other paths are derived from this automatically.
# Example in 01_run_preprocessing.R, 02_..., and 03_...
BASE_DIR <- "/path/to/your/project/on/hpc"

Data Structure: Make sure your raw data and configuration files are located in the subdirectories as defined in the SITES_CONFIG list in 01_run_preprocessing.R (e.g., 1_ilirney Lake/, A_aquatic_terrestrial/, etc.).
Execution on a SLURM Cluster
The included run_hpc.sh script is designed for the SLURM job scheduler. The execution is performed in three steps.
Start Stage 1 (Pre-processing):
Create a directory for log files and submit the job array. This command will start 6 jobs (0-5), one for each core.
mkdir -p logs
sbatch --array=0-5 run_hpc.sh stage1

Take note of the Job ID (e.g., 12345).
Start Stage 2 (Combining):
Submit this job after all jobs from Stage 1 have completed successfully. It will only start after the previous job array is finished.
# Replace <stage1_job_id> with the actual ID from the first submission
sbatch --dependency=afterok:12345 run_hpc.sh stage2

Take note of the new Job ID (e.g., 12346).
Start Stage 3 (Biomass Modeling):
Finally, submit the biomass modeling stage, which depends on the successful completion of the combination step.
# Replace <stage2_job_id> with the actual ID from the second submission
sbatch --dependency=afterok:12346 --array=0-5 run_hpc.sh stage3

After these steps are complete, the final biomass summary files will be located in the output directory specified by PATH_OUTPUT_DIR in the R scripts.


