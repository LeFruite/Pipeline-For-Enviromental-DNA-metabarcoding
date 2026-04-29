#!/bin/bash
set -e
set -u

# --- DEFAULT VARIABLES & CONFIGURATION ---
INTERACTIVE=true
START_STEP=1
TARGET_PATH=""
THREADS=45
DB_DIR="~/Documents/BioinformaticsVova/databases"

# Default Filtering & BLAST Thresholds
MIN_LEN=190
MAX_LEN=230
BLAST_MODE="local"
THRESH_PID=95
THRESH_EVAL=1e-6
THRESH_BITS=180

# --- HELP MESSAGE ---
usage() {
    echo "Usage: $0 [-s step_num] [-b] [input_file_or_directory]"
    echo "  -s  Start step (default: 1)"
    echo "  -b  Batch/Auto mode (Disables all interactive prompts, uses defaults)"
    echo "Steps:"
    echo "      1 = Quality Control (fastp)"
    echo "      2 = Conversion & Filter"
    echo "      3 = BLAST Alignment"
    echo "      4 = Taxonomy Mapping"
    echo "      5 = Final Filtering"
    exit 1
}

# --- ARGUMENT PARSING ---
while getopts "s:bh" opt; do
  case $opt in
    s) START_STEP=$OPTARG ;;
    b) INTERACTIVE=false ;; # Batch mode turns off interactivity
    h|*) usage ;;
  esac
done
shift $((OPTIND-1))

TARGET_PATH="${1:-}"

# --- INITIAL SETUP & PROMPTS (Only runs if INTERACTIVE=true) ---
echo "=========================================="
echo "      BIOINFORMATICS PIPELINE SETUP       "
echo "=========================================="

if [ -z "$TARGET_PATH" ]; then
    read -p "Enter the path to your input file OR directory: " TARGET_PATH
    if [ -z "$TARGET_PATH" ]; then
        echo "Error: No input provided."
        exit 1
    fi
fi

if [ "$INTERACTIVE" = true ]; then
    echo "Interactive mode enabled. Enter values or press [Enter] for defaults."
    read -p "Enter start step (1-5) [Default $START_STEP]: " input_step
    START_STEP=${input_step:-$START_STEP}

    read -p "Enter number of threads to use [Default $THREADS]: " input_threads
    THREADS=${input_threads:-$THREADS}

    if [ "$START_STEP" -le 4 ]; then
        read -p "Enter path to the database directory [Default $DB_DIR]: " input_db
        DB_DIR=${input_db:-$DB_DIR}
    fi

    if [ "$START_STEP" -le 2 ]; then
        read -p "Enter MINIMUM sequence length [Default $MIN_LEN]: " input_min
        MIN_LEN=${input_min:-$MIN_LEN}
        read -p "Enter MAXIMUM sequence length [Default $MAX_LEN]: " input_max
        MAX_LEN=${input_max:-$MAX_LEN}
    fi

    if [ "$START_STEP" -le 3 ]; then
        read -p "Use 'local' database or 'remote' NCBI? (local/remote) [Default $BLAST_MODE]: " input_blast
        BLAST_MODE=${input_blast:-$BLAST_MODE}
    fi

    if [ "$START_STEP" -le 5 ]; then
        read -p "Enter minimum Percent Identity [Default $THRESH_PID]: " input_pid
        THRESH_PID=${input_pid:-$THRESH_PID}
        read -p "Enter maximum E-value [Default $THRESH_EVAL]: " input_eval
        THRESH_EVAL=${input_eval:-$THRESH_EVAL}
        read -p "Enter minimum Bitscore [Default $THRESH_BITS]: " input_bits
        THRESH_BITS=${input_bits:-$THRESH_BITS}
    fi
else
    echo "Batch Mode (-b) Active: Using default parameters."
fi

# Set dynamic DB variables based on final DB_DIR
BLAST_DB="${DB_DIR}/bold_db"
TAXONOMY_MAP="${DB_DIR}/fasta_map.tsv"

# --- HELPER FUNCTIONS ---
confirm() {
    local step_name=$1
    if [ "$INTERACTIVE" = true ]; then
        echo -e "\n>>> FINISHED: $step_name"
        read -p "Press [Enter] to continue or [Ctrl+C] to exit..."
    fi
}

check_input() {
    local file_to_check=$1
    if [ ! -f "$file_to_check" ]; then
        echo "Error: Required input file '$file_to_check' does not exist! Skipping..."
        return 1
    fi
    return 0
}

# ==============================================================================
# MAIN PROCESSING FUNCTION (Runs for each file)
# ==============================================================================
process_sample() {
    local CURRENT_INPUT="$1"
    local BASENAME=$(basename "$CURRENT_INPUT")
    BASENAME="${BASENAME%.*}" 
    local RESULT_DIR="results_${BASENAME}"
    
    echo -e "\n=========================================="
    echo "PROCESSING SAMPLE: $BASENAME"
    echo "=========================================="
    
    mkdir -p "$RESULT_DIR"

    # Copy input file to results dir to preserve the original batch folder
    if [ -f "$CURRENT_INPUT" ]; then
        cp "$CURRENT_INPUT" "$RESULT_DIR/"
        CURRENT_INPUT="${RESULT_DIR}/$(basename "$CURRENT_INPUT")"
    fi

    # --- STEP 1: QUALITY CONTROL ---
    if [ "$START_STEP" -le 1 ]; then
        check_input "$CURRENT_INPUT" || return 1
        OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_clean.fastq"
        
        echo "[Step 1] Running fastp..."
        fastp -i "$CURRENT_INPUT" -o "$OUTPUT_FILE" -w "$THREADS" \
            --cut_front --cut_tail --cut_mean_quality 20 \
            --html "${RESULT_DIR}/${BASENAME}_fastp.html" --json "${RESULT_DIR}/${BASENAME}_fastp.json" > /dev/null 2>&1
        
        confirm "Quality Control"
        CURRENT_INPUT="$OUTPUT_FILE"
    fi

    # --- STEP 2: CONVERSION & LENGTH FILTERING ---
    if [ "$START_STEP" -le 2 ]; then
        check_input "$CURRENT_INPUT" || return 1
        TEMP_FASTA="${RESULT_DIR}/${BASENAME}_temp.fasta"
        OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_filtered.fasta"

        echo "[Step 2] Filtering sequences between ${MIN_LEN}bp and ${MAX_LEN}bp..."
        if [[ "$CURRENT_INPUT" == *.fastq || "$CURRENT_INPUT" == *.fq ]]; then
            seqtk seq -a "$CURRENT_INPUT" > "$TEMP_FASTA"
        else
            cp "$CURRENT_INPUT" "$TEMP_FASTA"
        fi
        
        seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" "$TEMP_FASTA" > "$OUTPUT_FILE"
        rm "$TEMP_FASTA"

        confirm "Seqkit Filtering"
        CURRENT_INPUT="$OUTPUT_FILE"
    fi

    # --- STEP 3: BLAST ALIGNMENT ---
    if [ "$START_STEP" -le 3 ]; then
        check_input "$CURRENT_INPUT" || return 1
        OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_blast.tsv"
        
        if [[ "$BLAST_MODE" == "remote" ]]; then
            echo "[Step 3] Running remote BLASTN against NCBI..."
            blastn -query "$CURRENT_INPUT" -db nt -remote -out "$OUTPUT_FILE" \
            -entrez_query '("CO1"[GENE] OR "COI"[GENE] OR "COX1"[GENE] OR "COXI"[GENE]) AND "Eukaryota"[ORGN] AND "BARCODE"[KYWD]' \
            -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdoms stitle" \
            -max_target_seqs 10
            
            CURRENT_INPUT="$OUTPUT_FILE"
            echo "Remote BLAST completed for $BASENAME. Stopping pipeline for this sample as requested."
            return 0 # Skips the rest of the steps for THIS file, moves to next file
        else
            echo "[Step 3] Running local BLASTN..."
            if [ ! -f "${BLAST_DB}.nal" ] && [ ! -f "${BLAST_DB}.nhr" ] && [ ! -f "${BLAST_DB}.nin" ]; then 
                echo "Error: Local BLAST DB not found at ${BLAST_DB}!"
                return 1
            fi
            
            blastn -query "$CURRENT_INPUT" -db "$BLAST_DB" -out "$OUTPUT_FILE" \
                -outfmt "6 qseqid sseqid pident length evalue bitscore" \
                -max_target_seqs 10 -num_threads "$THREADS"
        fi

        confirm "BLAST Alignment"
        CURRENT_INPUT="$OUTPUT_FILE"
    fi

    # --- STEP 4: TAXONOMY MAPPING ---
    if [ "$START_STEP" -le 4 ]; then
        check_input "$CURRENT_INPUT" || return 1
        OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_with_tax.tsv"
        
        echo "[Step 4] Mapping Taxonomy..."
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

    # --- STEP 5: FINAL FILTERING ---
    if [ "$START_STEP" -le 5 ]; then
        check_input "$CURRENT_INPUT" || return 1
        OUTPUT_FILE="${RESULT_DIR}/${BASENAME}_final.tsv"
        CONFLICTS_FILE="${RESULT_DIR}/${BASENAME}_conflicts.tsv"
        
        echo "[Step 5] Applying final quality filters (PID:${THRESH_PID}, EVAL:${THRESH_EVAL}, BITS:${THRESH_BITS})..."
        > "$OUTPUT_FILE"
        > "$CONFLICTS_FILE"

        awk -v p="$THRESH_PID" -v e="$THRESH_EVAL" -v b="$THRESH_BITS" \
            -v final_out="$OUTPUT_FILE" -v conflict_out="$CONFLICTS_FILE" \
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
            }' "$CURRENT_INPUT"

        confirm "Final Filtering"
    fi
    echo "Sample $BASENAME completed successfully."
}

# ==============================================================================
# EXECUTION LOGIC (Folder vs Single File)
# ==============================================================================

if [ -d "$TARGET_PATH" ]; then
    echo "Directory detected. Processing all files in $TARGET_PATH..."
    # Loops through any standard fastq/fasta extensions
    for file in "$TARGET_PATH"/*.{fastq,fq,fasta,fa}; do
        # Check if the glob actually found files
        if [ -f "$file" ]; then
            process_sample "$file"
        fi
    done
elif [ -f "$TARGET_PATH" ]; then
    echo "Single file detected."
    process_sample "$TARGET_PATH"
else
    echo "Error: Target '$TARGET_PATH' is neither a valid file nor a directory."
    exit 1
fi

echo -e "\n=== BATCH PIPELINE COMPLETE ==="
