# Pipeline-for-Environmental-DNA-(eDNA)-Metabarcoding
This pipeline automates data validation and taxonomic identification for sequencing data. The system is built using Bash, leveraging industry-standard bioinformatics tools for efficiency and speed.\
**You should use "pipeline.sh" file for completed product, down below will be command to install it**
## Overview
This is a step-by-step system that transforms raw FASTQ files into filtered files with completed taxonomic identifications.

### Step 1
The pipeline processes the raw FASTQ file and filters reads based on quality. While the Phred score threshold is customizable, it defaults to 20. This step also generates a summary statistics file for quality control review.
### Step 2
Converts FASTQ to FASTA format. It also filters sequences by length to remove outliers (reads that are too short or too long). The desired length range is user-defined.
### Step 3
Using BLASTn, the pipeline aligns sequences against a reference database. This supports both local and remote database queries.
### Step 4
The pipeline parses the .tsv alignment results to map specific taxonomy to the sequence data.
### Step 5
To ensure accuracy, the final results are filtered based on high-confidence metrics: % identity, e-value, length, and bit-score.

## Installation
You can download the pipeline script directly using curl:\
```curl -L https://raw.githubusercontent.com/LeFruite/Pipeline-For-Enviromental-DNA-metabarcoding/main/pipeline.sh -o pipeline.sh```

## Dependencies
Ensure you have the following tools installed on your system:\
```sudo apt update``` 
```sudo apt install fastp seqtk seqkit ncbi-blast+ r-base gawk coreutils```


## How to run this pipeline?
Before running the pipeline for the first time, you must grant the script execution permissions:\
```chmod +x pipeline.sh```
To execute the pipeline, run the following command followed by your input file:\
```./pipeline.sh Your_file.fastq```
