# Pipeline-For-Enviromental-DNA-metabarcoding
I have been developing a pipeline for data validation and taxonomic identification of sequencing data. This pipeline is utilising algorithms written in the Bash shell language.\
**You should use pipeline.sh for completed product, down below will be command to install it**
## What does it do exactly?
This is a step by step system that transforms the raw FastQ file into a file with completed taxonomic identification.

### Step 1
It takes FastaQ file and filters it by quality, you may change you desired Phred score, but it filters everything lower than 20 by default. During this step a file with useful statistics is generated, which you may look up.

### Step 2
It converts FastaQ into Fasta file. Also it filters by sequence length, to cut down short a long sequences. The length can be inputed.

### Step 3
Using BlastN we allign sequences in the file to the reference database. You can use both local and remote database.

### Step 4
Using .tsv of taxanomic allignment results, we take taxonomy and apply it to our results

### Step 5
It filters those results, so that only most probable are left. It filters by: % identity, e-value, length and score. 

## How to install this file?
Use this command to get my pipeline file.\
```curl -L https://raw.githubusercontent.com/LeFruite/Pipeline-For-Enviromental-DNA-metabarcoding/main/pipeline.sh -o pipeline.sh```
>or you can clone the whole folder

## How to run this pipeline?
Use this command to run my pipeline.\
```./pipeline Your_file.fastq```
