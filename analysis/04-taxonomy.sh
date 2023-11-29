#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script assigns taxonomy the de novo assembled genomes

cd $folder_shell
echo "04-taxonomy"

for i in {1..19}
do
    echo "$folder_raw/$batch_ids[$i]/$sample_ids[$i]"
    raw_reads="$folder_raw/$batch_ids[$i]/$sample_ids[$i]/reads/raw_reads.fastq.gz"
    filtered_reads="$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.fastq.gz"
    downsampled_reads="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/downsampled_reads.fastq"
    assembly_prefix="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/miniasm/draft_assembly"
    draft_genome="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/miniasm/draft_genome.fasta"
    downsampled_reads2="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/downsampled_reads2.fastq"

done

