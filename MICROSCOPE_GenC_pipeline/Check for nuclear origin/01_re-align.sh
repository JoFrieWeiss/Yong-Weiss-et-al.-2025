#!/bin/bash

#===========================================================================
# slurm batch script to taxon-specific re-mapping
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
#SBATCH --job-name=reM
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --array=1-9%9
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Josefine-Friederike.Weiss@awi.de

# ============ MODULES ============
# Load modules
module load bowtie2/2.5.1
module load samtools/1.20

#set -euo pipefail

# ============ CONFIGURATION ============
# parent
work_dir="/albedo/scratch/user/joweis001/BTOKO"
# path to out.metaDMG.subset
path2dmg="${work_dir}/output-lc/out.metaDMG.subset"
mapping_file="${work_dir}/sample2org2ref.tsv"
ref_base_dir="${work_dir}/references/ncbi_dataset/data"

# ============ GET SAMPLE ============
cd ${path2dmg}
sample_id=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))

# path to merge fastq.gz
fastq_end="_fastp_dedupe_merged.fq.gz"
fastq_input="${work_dir}/${sample_id}${fastq_end}"
output_dir="${work_dir}/alignment_output"

echo "========================================="
echo "Processing sample: $sample_id"
echo "========================================="

mkdir -p "${output_dir}/${sample_id}"

# ============ PROCESS EACH ORGANISM ============
# Read mapping file and process entries for this sample
while IFS=$'\t' read -r sample org ref_genome ref_size; do
  # Skip header
  if [[ "$sample" == "sample" ]]; then continue; fi
  
  # Only process current sample
  if [[ "$sample" != "$sample_id" ]]; then continue; fi
  
  echo ""
  echo "--- Processing: $org ($ref_genome) ---"
  
  # ============ STEP 1: Find seqid file ============
  seqid_file=$(find "${path2dmg}/${sample_id}" -name "*_${org}_*_seqid.txt" | head -1)
  
  if [ ! -f "$seqid_file" ]; then
    echo "ERROR: Seqid file not found for $org"
    continue
  fi
  
  echo "Step 1: Found seqid file: $(basename $seqid_file)"
  n_reads=$(wc -l < "$seqid_file")
  echo "  → $n_reads sequence IDs"
  
  # ============ STEP 2: Extract reads from fastq.gz ============
  fastq_subset="${output_dir}/${sample_id}/${sample_id}_${org}.fastq"
  
  if [ ! -f "$fastq_input" ]; then
    echo "ERROR: Input fastq not found: $fastq_input"
    continue
  fi
  
  echo "Step 2: Extracting reads from FASTQ..."
  
  # Method: grep with fixed strings (faster for many IDs)
  zcat "$fastq_input" | \
    awk 'NR==FNR {ids[$1]=1; next} 
         /^@/ {
           id=substr($1,2); 
           split(id,a,":"); 
           short_id=a[1]":"a[2]":"a[3]":"a[4]":"a[5]":"a[6]":"a[7];
           print_it = (short_id in ids)
         } 
         print_it' \
    "$seqid_file" - > "$fastq_subset"
  
  extracted=$(( $(wc -l < "$fastq_subset") / 4 ))
  echo "  → Extracted $extracted reads"
  
  if [ "$extracted" -eq 0 ]; then
    echo "WARNING: No reads extracted, skipping..."
    rm "$fastq_subset"
    continue
  fi
  
  # ============ STEP 3.1: Locate reference fasta  ============
  echo "Step 3.1: Locating reference for $ref_genome..."

  # Find gzipped reference
  #ref_fasta_original=$(find "$ref_base_dir" -maxdepth 1 -type f -name "${ref_genome}*.fna.gz" | head -1)
  ref_fasta_original=$(find "$ref_base_dir" -type f -name "${ref_genome}*.fna.gz" | head -1)

  if [ ! -f "$ref_fasta_original" ]; then
    # Try uncompressed
    #ref_fasta_original=$(find "$ref_base_dir" -maxdepth 1 -type f -name "${ref_genome}*.fna" | head -1)
    ref_fasta_original=$(find "$ref_base_dir" -type f -name "${ref_genome}*.fna" | head -1)
  fi

  if [ ! -f "$ref_fasta_original" ]; then
    echo "ERROR: Reference not found for: $ref_genome"
    echo "  Searched in: $ref_base_dir"
    rm "$fastq_subset"
    continue
  fi

  echo "  → Found: $(basename $ref_fasta_original)"

  # ============ STEP 3.2: Decompress if needed ============
  # Copy reference into sample-specific directory
  sample_ref_dir="${ref_base_dir}/sample_refs/${sample_id}"
  mkdir -p "$sample_ref_dir"
  ref_copy="${sample_ref_dir}/$(basename $ref_fasta_original)"

  if [ ! -f "$ref_copy" ]; then
    echo "Step 3.2: Copying reference to sample directory..."
    cp "$ref_fasta_original" "$ref_copy"
  else
    echo "Step 3.2: Reference already copied"
  fi

  # Only unzip inside sample-specific directory
  if [[ "$ref_copy" == *.gz ]]; then
    ref_fasta="${ref_copy%.gz}"
    if [ ! -f "$ref_fasta" ]; then
      echo "  → Decompressing within sample directory..."
      gunzip -c "$ref_copy" > "$ref_fasta"
    else
      echo "  → Uncompressed reference already exists"
    fi
  else
    ref_fasta="$ref_copy"
  fi
  echo "  → Using reference: $(basename $ref_fasta)"
  # ============ STEP 3.3: Build bowtie2 index inside sample directory ============
  index_base="${ref_fasta%.fna}"
  if [ ! -f "${index_base}.1.bt2" ]; then
    echo "Step 3.3: Building bowtie2 index..."
    srun bowtie2-build --threads 8 "$ref_fasta" "$index_base" </dev/null
  else
    echo "Step 3.3: Index already exists"
  fi

  # ============ STEP 4: Align reads ============
  bam_output="${output_dir}/${sample_id}/${sample_id}_${org}.bam"
  log_file="${output_dir}/${sample_id}/${sample_id}_${org}_bowtie2.log"
  
  echo "Step 4: Aligning reads with bowtie2..."
  
  srun bowtie2 \
    -x "$index_base" \
    -U "$fastq_subset" \
    -p 8 \
    --no-unal \
    </dev/null 2> "$log_file" | \
  samtools view -b - | \
  samtools sort -@ 4 -o "$bam_output"
  
  samtools index "$bam_output"
  
  # Get mapping stats
  mapped=$(samtools view -c -F 4 "$bam_output")
  total=$(samtools view -c "$bam_output")
  pct=$(awk -v m=$mapped -v t=$total 'BEGIN {printf "%.2f", (m/t)*100}')
  
  echo "  → Mapped: $mapped / $total ($pct%)"
  
  # ============ STEP 5: Cleanup ============
  echo "Step 5: Cleaning up intermediate files..."
  rm "$fastq_subset"
  
  echo "✓ Completed: $org"
  echo "  Output: $(basename $bam_output)"
  
  # ============ STEP 6: Generate Alignment Summary ============
  echo "Step 6: Generating alignment summary..."
  bam_output="${output_dir}/${sample_id}/${sample_id}_${org}.bam"
  log_file="${output_dir}/${sample_id}/${sample_id}_${org}_bowtie2.log"
  summary_file="${output_dir}/${sample_id}/${sample_id}_${org}_summary.txt"
  tsv_file="${output_dir}/${sample_id}/${sample_id}_${org}_alignments_counts.tsv"

  # Get basic stats
  total_alignments=$(samtools view -c "$bam_output")
  unique_reads=$(samtools view "$bam_output" | cut -f1 | sort -u | wc -l)

  # Save as TSV file (easy to import to Excel/R)
  {
    echo -e "reference_id\talignments\tsample\torganism\tref_genome"
    samtools view "$bam_output" | cut -f3 | sort | uniq -c | sort -rn | \
        awk -v sample="$sample_id" -v org="$org" -v ref="$ref_genome" \
        '{printf "%s\t%d\t%s\t%s\t%s\n", $2, $1, sample, org, ref}'
  } > "$tsv_file"

    echo "  → TSV saved to: $(basename $tsv_file)"

  echo "  → Total alignments: $total_alignments"
  echo "  → Unique reads: $unique_reads"

  # Count alignments per reference (accession.version)
  echo ""
  echo "  Top 20 aligned references:"
  samtools view "$bam_output" | cut -f3 | sort | uniq -c | sort -rn | head -20 | \
    awk '{printf "    %-40s %d alignments\n", $2, $1}'

  # Save summary to file
  {
    echo "========================================="
    echo "Alignment Summary for ${sample_id} - ${org}"
    echo "========================================="
    echo "Reference genome: $ref_genome"
    echo "Date: $(date)"
    echo ""
    echo "OVERALL STATISTICS:"
    echo "  Total alignments: $total_alignments"
    echo "  Unique reads: $unique_reads"
    echo "  Avg alignments per read: $(awk -v t=$total_alignments -v u=$unique_reads 'BEGIN {printf "%.2f", t/u}')"
    echo ""
    echo "ALIGNMENTS BY REFERENCE ACCESSION:"
    samtools view "$bam_output" | cut -f3 | sort | uniq -c | sort -rn | awk '{printf "  %-30s %10d alignments\n", $2, $1}'
    echo ""
    echo "TOP 10 MOST ALIGNED REFERENCES:"
    samtools view "$bam_output" | cut -f3 | sort | uniq -c | sort -rn | head -10 | awk '{printf "  %2d. %-30s %10d alignments\n", NR, $2, $1}'
  } < /dev/null > "$summary_file"

    echo "  → Summary saved to: $(basename $summary_file)"

done < <(tail -n +2 "$mapping_file")

echo ""
echo "========================================="
echo "Sample $sample_id complete!"
echo "========================================="
echo "Output directory: ${output_dir}/${sample_id}"
ls -lh "${output_dir}/${sample_id}"
