#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

folder_ibd=$folder_genomics/ibd # ibd folder
folder_refs=$folder_ibd/refs # reference folder
folder_query=$folder_ibd/query # query folder. only do S meliloti and medicae
folder_blast=$folder_ibd/blast_results
mkdir -p $folder_refs $folder_blast $folder_query

# Download ref genomes
mamba activate ncbi-datasets
# S. meliloti
datasets download genome accession GCF_000006965.1 --include genome \
  --filename $folder_refs/s_meliloti.zip
# S. medicae
datasets download genome accession GCF_000017145.1 --include genome \
  --filename $folder_refs/s_medicae.zip

# Extract FASTA files
unzip -o $folder_refs/s_meliloti.zip -d $folder_refs/s_meliloti
unzip -o $folder_refs/s_medicae.zip -d $folder_refs/s_medicae

meliloti_fna=$(find "$folder_refs/s_meliloti" -name "*.fna" | head -n1)
medicae_fna=$(find "$folder_refs/s_medicae" -name "*.fna" | head -n1)

# Build on single blast database
mamba activate blast
cat $meliloti_fna $medicae_fna > $folder_refs/Sinorhizobium_combined_refs.fna
makeblastdb -in $folder_refs/Sinorhizobium_combined_refs.fna -dbtype nucl

# Move query fasta
for i in 4 5 6 8 9 10 11 13 16 17 19 20 {21..27} {29..37} {39..45};
do
    cp $folder_genomics/fasta/genomes/g$i.fasta $folder_query/g$i.fasta
done

# Blast
for query in $folder_query/*.fasta(.N); do
    sample=$(basename "$query" .fasta)
    echo "  - Processing $sample"

    blastn \
      -query $query\
      -db $folder_refs/Sinorhizobium_combined_refs.fna \
      -out $folder_blast/${sample}_vs_refs.tsv \
      -outfmt "6 qseqid sseqid qlen pident length evalue bitscore" \
      -max_target_seqs 3 \
      -num_threads 8
done
