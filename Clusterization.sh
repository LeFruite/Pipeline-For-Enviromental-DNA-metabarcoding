python -m venv speconsense-env
source speconsense-env/bin/activate  # On Windows: speconsense-env\Scripts\activate

# --min-size 10 --algorithm greedy --min-identity 0.85
speconsense input.fastq -O clust_results/ --threads 8

# Merge clusters with >= 99% identity (1% dissimilarity)
speconsense-summarize --source . --summary-dir MyResults --min-identity 0.99

# Install via Conda
conda install -c bioconda taxonkit

# Download and set up the NCBI taxonomy data (approx 400MB)
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
mkdir -p ~/.taxonkit
mv names.dmp nodes.dmp delnodes.dmp merged.dmp ~/.taxonkit

taxonkit lineage -i 7 blast_results.tsv | taxonkit reformat -i 11 -f "{k};{p};{c};{o};{f};{g};{s}" > blast_with_taxonomy.tsv
