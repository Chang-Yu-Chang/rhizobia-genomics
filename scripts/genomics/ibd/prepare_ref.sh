#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

folder_ibd=$folder_genomics/ibd
folder_refs=$folder_ibd/refs
mkdir -p $folder_refs $folder_ibd/blast_results


# Download ref genomes
mamba activate ncbi-datasets
# S. meliloti
datasets download genome taxon "Sinorhizobium meliloti" \
  --reference --include genome \
  --filename $folder_refs/s_meliloti.zip

# S. medicae
datasets download genome taxon "Sinorhizobium medicae" \
  --reference --include genome \
  --filename $folder_refs/s_medicae.zip

# Extracting FASTA files
unzip -o $folder_refs/s_meliloti.zip -d $folder_refs/s_meliloti
unzip -o $folder_refs/s_medicae.zip -d $folder_refs/s_medicae

meliloti_fna=$(find "$folder_refs/s_meliloti" -name "*.fna" | head -n1)
medicae_fna=$(find "$folder_refs/s_medicae" -name "*.fna" | head -n1)

cp $meliloti_fna $folder_refs/Sinorhizobium_meliloti_ref.fna
cp $medicae_fna $folder_refs/Sinorhizobium_medicae_ref.fna


# ========= 3. BUILD A SINGLE BLAST DATABASE ==========
mamba activate blast
cat "$folder_refs/Sinorhizobium_meliloti_ref.fna" \
    "$folder_refs/Sinorhizobium_medicae_ref.fna" \
    > "$folder_refs/Sinorhizobium_combined_refs.fna"
makeblastdb -in "$folder_refs/Sinorhizobium_combined_refs.fna" -dbtype nucl


for query in "$folder_query"/*.fasta(.N); do
    sample=$(basename "$query" .fasta)
    echo "  - Processing $sample"

    blastn \
      -query "$query" \
      -db "$folder_refs/Sinorhizobium_combined_refs.fna" \
      -out "$folder_blast/${sample}_vs_combined.tsv" \
      -outfmt "6 qseqid sseqid pident length evalue bitscore" \
      -max_target_seqs 3 \
      -num_threads 8
done
