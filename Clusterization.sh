python -m venv speconsense-env
source speconsense-env/bin/activate  # On Windows: speconsense-env\Scripts\activate

# --min-size 10 --algorithm greedy --min-identity 0.85
speconsense input.fastq -O clust_results/ --threads 8

# Merge clusters with >= 99% identity (1% dissimilarity)
speconsense-summarize --source . --summary-dir MyResults --min-identity 0.99
