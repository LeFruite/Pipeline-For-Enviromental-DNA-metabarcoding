#!/bin/bash
set -e
set -u

# --- CONFIGURATION ---
# Сurrently this part of code is custom made for my PC. If you want to use it, you probably need to configure it yourself
DB_DIR="./databases" # This and next 2 variables are required to do allignment and taxonomic identification, if you are using local database
BLAST_DB="${DB_DIR}/bold_db"
TAXONOMY_MAP="${DB_DIR}/fasta_map.tsv"
#R_SCRIPT="./visualize.R" #This code is for visualisation after the taxanomic identification was made. WIP 
THREADS=45 #Number of threads, speeds up allignment. Need to make it customisable

# --- DEFAULT VARIABLES ---
INTERACTIVE=true #Changed it to be always true, might change it back
START_STEP=""
CURRENT_INPUT=""

# --- HELP MESSAGE ---
usage() {
    echo "Usage: $0 [-s step_num] [input_file]"
    echo "  -s  Start step (default: 1):"
    echo "      1 = Quality Control (fastp)      [Input: FASTQ]"
    echo "      2 = Conversion & Filter          [Input: Clean FASTQ]"
    echo "      3 = BLAST Alignment              [Input: FASTA]"
    echo "      4 = Taxonomy Mapping             [Input: BLAST TSV]"
    echo "      5 = Final Filtering              [Input: Mapped TSV]"
    echo "      6 = Visualization (R)            [Input: Final TSV]"
    echo "Note: Script is fully interactive and will prompt for missing inputs."
    exit 1
}

# --- ARGUMENT PARSING ---
while getopts "s:h" opt; do
  case $opt in
    s) START_STEP=$OPTARG ;;
    h|*) usage ;;
  esac
done
shift $((OPTIND-1))

# --- INTERACTIVE PROMPTS FOR MISSING ARGS ---
echo "=========================================="
echo "      BIOINFORMATICS PIPELINE SETUP       "
echo "=========================================="

if [ -z "$START_STEP" ]; then
    read -p "Enter start step (1-6) [Default 1]: " START_STEP
    START_STEP=${START_STEP:-1}
fi

if [ -n "${1:-}" ]; then
    CURRENT_INPUT="$1"
else
    read -p "Enter the path to your input file: " CURRENT_INPUT
    if [ -z "$CURRENT_INPUT" ]; then
        echo "Error: No input file provided."
        exit 1
    fi
fi

# --- DIRECTORY SETUP ---
BASENAME=$(basename "$CURRENT_INPUT")
BASENAME="${BASENAME%.*}" 
RESULT_DIR="results_${BASENAME}"
mkdir -p "$RESULT_DIR"

echo "=========================================="
echo "Start Step: $START_STEP"
echo "Input File: $CURRENT_INPUT"
echo "Output Dir: $RESULT_DIR"
echo "=========================================="

# --- HELPER FUNCTIONS ---
confirm() {
    local step_name=$1
    if [ "$INTERACTIVE" = true ]; then
        echo -e "\n>>> FINISHED: $step_name"
        read -p "Press [Enter] to continue to the next step or [Ctrl+C] to exit..."
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
    CURRENT_INPUT="$OUTPUT_FILE"
fi

# ==============================================================================
# STEP 2: CONVERSION & LENGTH FILTERING
# ==============================================================================
if [ "$START_STEP" -le 2 ]; then
    check_input
    TEMP_FASTA="${RESULT_DIR}/${BASENAME}_temp.fasta"
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_filtered.fasta"

    echo -e "\n[Step 2/6] Sequence Conversion and Filtering"
    
    # Interactive Length Inputs
    read -p "Enter MINIMUM sequence length [Default 600]: " MIN_LEN
    MIN_LEN=${MIN_LEN:-600}
    
    read -p "Enter MAXIMUM sequence length [Default 720]: " MAX_LEN
    MAX_LEN=${MAX_LEN:-720}

    echo "Filtering sequences between ${MIN_LEN}bp and ${MAX_LEN}bp..."

    if [[ "$CURRENT_INPUT" == *.fastq ]]; then
        seqtk seq -a "$CURRENT_INPUT" > "$TEMP_FASTA"
    else
        cp "$CURRENT_INPUT" "$TEMP_FASTA"
    fi
    
    seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" "$TEMP_FASTA" > "$OUTPUT_FILE"
    rm "$TEMP_FASTA"

    confirm "Seqkit Filtering"
    CURRENT_INPUT="$OUTPUT_FILE"
fi

# ==============================================================================
# STEP 3: BLAST ALIGNMENT
# I am still coding a way to make taxonomy for NCBI database, so if you choose remote
# analysis program will end your analysis after this step. 
# ==============================================================================
if [ "$START_STEP" -le 3 ]; then
    check_input
    OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_blast.tsv"
    
    echo -e "\n[Step 3/6] BLAST Alignment Setup"
    
    # Interactive BLAST Target Setup
    read -p "Use 'local' database or 'remote' NCBI? (local/remote) [Default local]: " BLAST_MODE
    BLAST_MODE=${BLAST_MODE:-local}
    
    if [[ "$BLAST_MODE" == "remote" ]]; then
        echo "Running remote BLASTN against NCBI 'nt' database (This may take a while)..."
        blastn -query "$CURRENT_INPUT" -db nt -remote -out "$OUTPUT_FILE" \
            -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdoms stitle" \
            -max_target_seqs 10
    else
        echo "Running local BLASTN..."
        if [ ! -f "${BLAST_DB}.nal" ] && [ ! -f "${BLAST_DB}.nhr" ] && [ ! -f "${BLAST_DB}.nin" ]; then 
            echo "Error: Local BLAST DB not found at ${BLAST_DB}!"
            exit 1
        fi
        
        blastn -query "$CURRENT_INPUT" -db "$BLAST_DB" -out "$OUTPUT_FILE" \
            -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdoms stitle" \
            -max_target_seqs 10 -num_threads "$THREADS"
    fi

    confirm "BLAST Alignment"
    CURRENT_INPUT="$OUTPUT_FILE"

    # Termination Logic for Remote BLAST
    if [[ "$BLAST_MODE" == "remote" ]]; then
        echo -e "\n=========================================="
        echo "Remote BLAST completed."
        echo "Results saved in: $CURRENT_INPUT"
        echo "Terminating pipeline as requested for remote runs."
        echo "=========================================="
        exit 0
    fi
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
    
    echo -e "\n[Step 5/6] Final Quality Filtering"
    
    # Interactive Filtering Thresholds
    read -p "Enter minimum Percent Identity (Field 3) [Default 90]: " THRESH_PID
    THRESH_PID=${THRESH_PID:-90}

    read -p "Enter maximum E-value (Field 5) [Default 1e-18]: " THRESH_EVAL
    THRESH_EVAL=${THRESH_EVAL:-1e-18}

    read -p "Enter minimum Bitscore (Field 6) [Default 120]: " THRESH_BITS
    THRESH_BITS=${THRESH_BITS:-120}

    read -p "Enter minimum threshold for Field 10 [Default 80]: " THRESH_F10
    THRESH_F10=${THRESH_F10:-80}

    echo "Applying filters: Identity >= ${THRESH_PID}, E-value <= ${THRESH_EVAL}, Bitscore > ${THRESH_BITS}, Field 10 >= ${THRESH_F10}..."

    # Deduplicate 
    awk -F'\t' '!seen[$1]++' "$CURRENT_INPUT" > "$TEMP_UNIQUE"
    
    # Quality Filter passing user variables to awk
    awk -v p="$THRESH_PID" -v e="$THRESH_EVAL" -v b="$THRESH_BITS" -v f10="$THRESH_F10" \
        -F'\t' '($3 >= p) && ($5 <= e) && ($6 > b) && ($10 >= f10)' "$TEMP_UNIQUE" > "$OUTPUT_FILE"

    confirm "Final Filtering"
    CURRENT_INPUT="$OUTPUT_FILE"
fi

# ==============================================================================
	
# STEP 6: VISUALIZATION (R)
	
# ==============================================================================
#	
#if [ "$START_STEP" -le 6 ]; then
#	
#    check_input
#	
#    echo -e "\n[Step 6/6] Running R Visualization..."
#	
#
#	
#    if command -v Rscript &> /dev/null; then
#	
#        Rscript "$R_SCRIPT" "$CURRENT_INPUT" "$RESULT_DIR"
#	
#    else
#	
#        echo "Warning: Rscript not found, skipping visualization."
#	
#    fi
#	
#fi
echo -e "\n=== PIPELINE COMPLETE ==="
echo "Final results in: $RESULT_DIR"
