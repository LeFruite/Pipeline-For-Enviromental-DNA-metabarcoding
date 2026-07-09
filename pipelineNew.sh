#!/bin/bash
set -e
set -u

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Define all variables here at the top for easy editing
wd="./"
input_file=${1:-}

# Hardware and Databases
threads=4
db_dir="./databases"
blast_db="${db_dir}/bold_db"
taxonomy_map="${db_dir}/fasta_map.tsv"

# Filtering Thresholds
min_len=600
max_len=720
thresh_pid=95
thresh_eval=1e-6
thresh_bits=180

# ==============================================================================
# INITIALIZATION
# ==============================================================================
cd "$wd"

if [ -z "$input_file" ]; then
    echo "Error: No input file provided."
    echo "Usage: $0 <input_file.fastq>"
    exit 1
fi

basename=$(basename "$input_file" | sed 's/\.[^.]*$//')
result_dir="results_${basename}"
mkdir -p "$result_dir"
mkdir -p logs

echo " "
echo "===================================="
echo "Processing sample: $basename"
echo "===================================="

# ==============================================================================
# PIPELINE EXECUTION
# ==============================================================================

# --- STEP 1: Quality Control ---
echo "-- quality control (fastp)"
clean_fastq="${result_dir}/${basename}_clean.fastq"

fastp -i "$input_file" -o "$clean_fastq" -w "$threads" \
    --cut_front --cut_tail --cut_mean_quality 20 \
    --html "${result_dir}/${basename}_fastp.html" \
    --json "${result_dir}/${basename}_fastp.json" 2> logs/fastp.log


# --- STEP 2: Conversion & Length Filtering ---
echo "-- sequence conversion and filtering"
temp_fasta="${result_dir}/${basename}_temp.fasta"
filtered_fasta="${result_dir}/${basename}_filtered.fasta"

if [[ "$clean_fastq" == *.fastq ]]; then
    seqtk seq -a "$clean_fastq" > "$temp_fasta"
else
    cp "$clean_fastq" "$temp_fasta"
fi

seqkit seq -m "$min_len" -M "$max_len" "$temp_fasta" > "$filtered_fasta"
rm "$temp_fasta"


# --- STEP 3: BLAST Alignment ---
echo "-- BLAST alignment (local)"
blast_tsv="${result_dir}/${basename}_blast.tsv"

if [ ! -f "${blast_db}.nhr" ] && [ ! -f "${blast_db}.nal" ]; then 
    echo "Error: Local BLAST DB not found at ${blast_db}!"
    exit 1
fi

blastn -query "$filtered_fasta" -db "$blast_db" -out "$blast_tsv" \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 10 -num_threads "$threads" 2> logs/blast.log


# --- STEP 4: Taxonomy Mapping ---
echo "-- mapping taxonomy"
mapped_tsv="${result_dir}/${basename}_with_tax.tsv"

awk 'BEGIN { FS="\t"; OFS="\t" }
NR==FNR {
    sub(/^>/, "", $1);
    split($1, parts, /[ \t]/);
    id = parts[1];
    val = $0; sub(/^[^ \t]+[ \t]+/, "", val);
    map[id] = val;
    next
}
{
    split($2, a, "|");
    id = a[1];
    tax = (id in map ? map[id] : "NA");
    print $0, tax
}' "$taxonomy_map" "$blast_tsv" > "$mapped_tsv"


# --- STEP 5: Final Filtering & Congruency Check ---
echo "-- final filtering & congruency check"
final_tsv="${result_dir}/${basename}_final.tsv"
conflicts_tsv="${result_dir}/${basename}_conflicts.tsv"

> "$final_tsv"
> "$conflicts_tsv"

awk -v p="$thresh_pid" -v e="$thresh_eval" -v b="$thresh_bits" \
    -v final_out="$final_tsv" -v conflict_out="$conflicts_tsv" \
    'BEGIN { FS="\t"; OFS="\t" }
    {
        if ($3 >= p && $5 <= e && $6 > b) {
            qseqid = $1
            if (!(qseqid in count)) {
                order[++num_queries] = qseqid
                count[qseqid] = 0
            }
            count[qseqid]++
            idx = count[qseqid]
            lines[qseqid, idx] = $0
            taxonomies[qseqid, idx] = $7 
        }
    }
    END {
        for (i = 1; i <= num_queries; i++) {
            qseqid = order[i]
            c = count[qseqid]

            if (c == 1) {
                print lines[qseqid, 1] > final_out
            } else {
                conflict = 0
                base_tax = taxonomies[qseqid, 1]

                for (j = 2; j <= c; j++) {
                    if (taxonomies[qseqid, j] != base_tax) {
                        conflict = 1
                        break
                    }
                }

                if (conflict == 0) {
                    print lines[qseqid, 1] > final_out
                } else {
                    for (j = 1; j <= c; j++) {
                        print lines[qseqid, j] > conflict_out
                    }
                }
            }
        }
    }' "$mapped_tsv"

echo " "
echo "===================================="
echo "Pipeline complete!"
echo "Results saved in: $result_dir"
echo "===================================="