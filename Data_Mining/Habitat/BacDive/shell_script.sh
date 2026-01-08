#!/bin/bash

#SBATCH --account=your_account_here
#SBATCH --job-name=database_scraping
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL        # Notify on start, end, and fail
#SBATCH --mail-user=your_email@example.com

python bacdive_scraping_v2.1.py
