#!/bin/bash
set -e
set -u

# --- CONFIGURATION ---
DB_DIR="./databases"
BLAST_DB="${DB_DIR}/bold_db"
TAXONOMY_MAP="${DB_DIR}/fasta_map.tsv"
R_SCRIPT="./visualize.R"
THREADS=45

# --- DEFAULT VARIABLES ---
INTERACTIVE=false
START_STEP=1
INPUT_FILE=""

# --- HELP MESSAGE ---
usage() {
    echo "Usage: $0 [-i] [-s step_num] <input_file>"
    echo "  -i  Interactive mode (ask before each step)"
    echo "  -s  Start step (default: 1):"
    echo "      1 = Quality Control (fastp)      [Input: FASTQ]"
    echo "      2 = Conversion & Filter          [Input: Clean FASTQ]"
    echo "      3 = BLAST Alignment              [Input: FASTA]"
    echo "      4 = Taxonomy Mapping             [Input: BLAST TSV]"
    echo "      5 = Final Filtering              [Input: Mapped TSV]"
    echo "      6 = Visualization (R)            [Input: Final TSV]"
    echo "Example: $0 -s 4 my_blast_results.tsv"
    exit 1
}

# --- ARGUMENT PARSING ---
while getopts "is:" opt; do
  case $opt in
    i) INTERACTIVE=true ;;
    s) START_STEP=$OPTARG ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${1:-}" ]; then
    echo "Error: No input file provided."
    usage
fi

# The file provided by the user is the starting point
CURRENT_INPUT="$1"

# --- DIRECTORY SETUP ---
# We derive the basename from the input file to name the results folder
# e.g., if input is "data/sample1.fastq", basename is "sample1"
# e.g., if input is "sample1_blast.tsv", basename is "sample1_blast"
BASENAME=$(basename "$CURRENT_INPUT")
BASENAME="${BASENAME%.*}" 
RESULT_DIR="results_${BASENAME}"
mkdir -p "$RESULT_DIR"

echo "=========================================="
echo "PIPELINE CONFIGURATION"
echo "Start Step: $START_STEP"
echo "Input File: $CURRENT_INPUT"
echo "Output Dir: $RESULT_DIR"
echo "=========================================="

# --- HELPER FUNCTIONS ---
confirm() {
    local step_name=$1
    if [ "$INTERACTIVE" = true ]; then
        echo -e "\n>>> FINISHED: $step_name"
        read -p "Press [Enter] to continue or [Ctrl+C] to exit..."
    fi
}

check_input() {
    if [ ! -f "$CURRENT_INPUT" ]; then
        echo "Error: Input file '$CURRENT_INPUT' for this step does not exist!"
        exit 1
    fi
}

# ==============================================================================
# STEP 1: QUALITY CONTROL (fastp)
# ==============================================================================
if [ "$START_STEP" -le 1 ]; then
    check_input
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_clean.fastq"
    
    echo -e "\n[Step 1/6] Running fastp..."
    fastp -i "$CURRENT_INPUT" -o "$OUTPUT_FILE" -w "$THREADS" \
        --cut_front --cut_tail --cut_mean_quality 20 \
        --html "${RESULT_DIR}/${BASENAME}_fastp.html" --json "${RESULT_DIR}/${BASENAME}_fastp.json" > /dev/null 2>&1
    
    confirm "Quality Control"
    CURRENT_INPUT="$OUTPUT_FILE" # Pass output to next step
fi

# ==============================================================================
# STEP 2: CONVERSION & LENGTH FILTERING
# ==============================================================================
if [ "$START_STEP" -le 2 ]; then
    check_input
    TEMP_FASTA="${RESULT_DIR}/${BASENAME}_temp.fasta"
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_filtered.fasta"

    echo -e "\n[Step 2/6] Converting and Filtering (600-720bp)..."
    # Detect if input is FASTQ or FASTA
    if [[ "$CURRENT_INPUT" == *.fastq ]]; then
        seqtk seq -a "$CURRENT_INPUT" > "$TEMP_FASTA"
    else
        cp "$CURRENT_INPUT" "$TEMP_FASTA"
    fi
    
    seqkit seq -m 600 -M 720 "$TEMP_FASTA" > "$OUTPUT_FILE"
    rm "$TEMP_FASTA"

    confirm "Seqkit Filtering"
    CURRENT_INPUT="$OUTPUT_FILE"
fi
# ==============================================================================
# STEP 3: BLAST ALIGNMENT
# ==============================================================================
if [ "$START_STEP" -le 3 ]; then
    check_input
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_blast.tsv"
    
    echo -e "\n[Step 3/6] Running BLASTN..."
    if [ ! -f "${BLAST_DB}.nal" ]; then echo "Error: BLAST DB not found!"; exit 1; fi
    
    blastn -query "$CURRENT_INPUT" -db "$BLAST_DB" -out "$OUTPUT_FILE" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdoms stitle" \
        -max_target_seqs 10 -num_threads "$THREADS"

    confirm "BLAST Alignment"
    CURRENT_INPUT="$OUTPUT_FILE"
fi

# ==============================================================================
# STEP 4: TAXONOMY MAPPING
# ==============================================================================
if [ "$START_STEP" -le 4 ]; then
    check_input
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_with_tax.tsv"
    
    echo -e "\n[Step 4/6] Mapping Taxonomy..."
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
    }' "$TAXONOMY_MAP" "$CURRENT_INPUT" > "$OUTPUT_FILE"

    confirm "Taxonomy Mapping"
    CURRENT_INPUT="$OUTPUT_FILE"
fi

# ==============================================================================
# STEP 5: FINAL FILTERING
# ==============================================================================
if [ "$START_STEP" -le 5 ]; then
    check_input
    TEMP_UNIQUE="${RESULT_DIR}/${BASENAME}_unique.tsv"
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_final.tsv"
    
    echo -e "\n[Step 5/6] Final Filtering..."
    # Deduplicate
    awk -F'\t' '!seen[$1]++' "$CURRENT_INPUT" > "$TEMP_UNIQUE"
    # Quality Filter
    awk -F'\t' '($3 >= 90) && ($5 <= 1e-18) && ($6 > 120) && ($10 >= 80)' "$TEMP_UNIQUE" > "$OUTPUT_FILE"

    confirm "Final Filtering"
    CURRENT_INPUT="$OUTPUT_FILE"
fi

# ==============================================================================
# STEP 6: VISUALIZATION (R)
# ==============================================================================
if [ "$START_STEP" -le 6 ]; then
    check_input
    echo -e "\n[Step 6/6] Running R Visualization..."
    
    if command -v Rscript &> /dev/null; then
        Rscript "$R_SCRIPT" "$CURRENT_INPUT" "$RESULT_DIR"
    else
        echo "Warning: Rscript not found, skipping visualization."
    fi
fi

echo -e "\n=== PIPELINE COMPLETE ==="
echo "Final results in: $RESULT_DIR"
