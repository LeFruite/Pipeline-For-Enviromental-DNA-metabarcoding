#!/bin/bash
set -e
set -u

# ==============================================================================
# CONFIGURATION
# ==============================================================================
wd="./"
input_fastq=${1:-}

# Hardware
threads=8

# Database Requirements
# WARNING: Your database FASTA headers MUST be in the SINTAX format:
# Example: >seq1;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales...
sintax_db="./databases/reference_sintax_db.fasta"

# Nanopore-Specific Parameters
# Lower the SINTAX cutoff (default is usually 0.8 for Illumina, 0.6 is safer for ONT)
sintax_cutoff=0.6 
min_len=1300      # Adjust based on your expected amplicon size (e.g., 1.5kb for 16S)
max_len=1700

# ==============================================================================
# INITIALIZATION
# ==============================================================================
cd "$wd"

if [ -z "$input_fastq" ]; then
    echo "Error: No input FASTQ file provided."
    echo "Usage: $0 <input_file.fastq>"
    exit 1
fi

basename=$(basename "$input_fastq" | sed 's/\.[^.]*$//')
result_dir="results_vsearch_${basename}"
mkdir -p "$result_dir"
mkdir -p logs

echo " "
echo "===================================="
echo "Processing Nanopore sample: $basename"
echo "===================================="

# --- STEP 1: Quality & Length Filtering ---
# Nanopore data often contains short junk reads or ultra-long chimeric reads.
# We filter them and convert FASTQ to FASTA simultaneously.
echo "-- Filtering length and converting to FASTA"
filtered_fasta="${result_dir}/${basename}_filtered.fasta"

vsearch --fastq_filter "$input_fastq" \
    --fastq_minlen "$min_len" \
    --fastq_maxlen "$max_len" \
    --fastaout "$filtered_fasta" \
    --threads "$threads" 2> logs/vsearch_filter.log

# --- STEP 2: SINTAX Taxonomic Assignment ---
echo "-- Running SINTAX classification"
sintax_out="${result_dir}/${basename}_sintax.tsv"

if [ ! -f "$sintax_db" ]; then
    echo "Error: SINTAX database not found at ${sintax_db}!"
    exit 1
fi

# We use '--strand both' because Nanopore reads can be in either orientation
vsearch --sintax "$filtered_fasta" \
    --db "$sintax_db" \
    --tabbedout "$sintax_out" \
    --sintax_cutoff "$sintax_cutoff" \
    --strand both \
    --threads "$threads" 2> logs/vsearch_sintax.log

# --- STEP 3: Clean up and Format Results ---
echo "-- Formatting taxonomy results"
final_taxonomy="${result_dir}/${basename}_taxonomy_clean.tsv"

# The standard SINTAX output contains 4 columns:
# 1. Query ID | 2. Predicted taxonomy | 3. Strand | 4. Taxonomy of the top hit
# This AWK command grabs the ID and the predicted taxonomy, stripping out SINTAX formatting characters.
awk -F'\t' '{
    if ($2 != "") {
        clean_tax = $2;
        gsub(/[a-z]:/, "", clean_tax); # Removes the d:, p:, c: prefixes
        gsub(/"/, "", clean_tax);      # Removes quotes
        print $1 "\t" clean_tax;
    } else {
        print $1 "\tUnclassified";
    }
}' "$sintax_out" > "$final_taxonomy"

echo " "
echo "===================================="
echo "VSEARCH SINTAX complete!"
echo "Raw output:   $sintax_out"
echo "Clean output: $final_taxonomy"
echo "===================================="