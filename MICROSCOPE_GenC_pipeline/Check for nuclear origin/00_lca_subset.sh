#!/bin/bash

#===========================================================================
# slurm batch script to extract taxon-specific reads from lca file
# on several sample using arrays
# Version 0.0 
#
# by Sisi Liu
# 
# contact: sisi.liu@awi.de
#
# slurm options and variables under >set required variables< 
# have to be modified by the user
#=============================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=exR
#SBATCH --partition=smp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --array=1-2%2
#SBATCH --mem=10G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Josefine-Friederike.Weiss@awi.de

# ============ CONFIG ============
# parent
work_dir="/albedo/scratch/user/joweis001/BTOKO"
# path to out.metaDMG
path2dmg="/albedo/work/projects/p_biodiv_shotgun/04_HOLI-post/genC"
# get sample list and past to slurm
cd ${path2dmg}
#sample_id=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))
sample_id="JK079L-9_S9"

# which min len
len_thd="L30"
# path to lca gz file
path2lca="${path2dmg}/${sample_id}/data/lca"
# lca gz
lca_gz="${path2lca}/${sample_id}_${len_thd}.lca.txt.gz"
# output
out_dir="${work_dir}/output-lc/out.metaDMG.subset/${sample_id}"
mkdir -p ${out_dir}

# processing
# zgrep '3688:Salicaceae:"family"' "${lca_gz}" | grep -v '^#' \
#   > ${out_dir}/JK081L-11_S11_L30_${prefix}.txt

# Define your patterns
patterns=(
  '3688:Salicaceae:"family"'
  '48230:Dryas:"genus"'
  '4345:Ericaceae:"family"'
  '8015:Salmonidae:"family"'
  '9376:Soricidae:"family"'
  '9976:Ochotonidae:"family"'
  '9655:Mustelidae:"family"'
)

echo "Extracting taxa from lca gz file..."

# Process each pattern
for pattern in "${patterns[@]}"; do
  # Remove quotes from pattern for parsing
  clean_pattern="${pattern//\"/}"
  
  # Parse pattern using IFS
  IFS=":" read -ra parts <<< "$clean_pattern"
  taxid="${parts[0]}"
  org="${parts[1]}"
  rank="${parts[2]}"
  
  # Create sequence ID for filename
  prefix="${taxid}_${org}_${rank}_lca"
  
  # Create output file
  sub_lca="${out_dir}/${sample_id}_${prefix}.txt"
  
  # subset lca
  zgrep "$pattern" "${lca_gz}" | grep -v '^#' > "$sub_lca"

  # get the seq id only
  prefix="${taxid}_${org}_${rank}_seqid"
  sub_lca_seq="${out_dir}/${sample_id}_${prefix}.txt"
  cut -f1 "$sub_lca" | cut -d':' -f1-7 > "$sub_lca_seq"

  # Report
  count=$(wc -l < "$sub_lca")
  echo "  ${org} (${rank} ${taxid}): $count reads → $(basename $sub_lca)"
  
  count=$(wc -l < "$sub_lca_seq")
  echo "  ${org} (${rank} ${taxid}): $count seqid → $(basename $sub_lca_seq)"
done

echo "step 0: lca subset, Done! Files created in: $out_dir"
