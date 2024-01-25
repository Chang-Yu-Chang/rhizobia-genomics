#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script performs quality control on the raw nanapore whole-genome long reads

cd $folder_shell

for i in {1..38}
do
    echo $genome_ids[$i]
    raw_reads="$folder_genomics/raw_reads/$genome_ids[$i].fastq.gz"
    filtered_reads="$folder_genomics/assembly/$genome_ids[$i]/filtered_reads.fastq.gz"

    # Remove the worst 5% reads via filtlong
    mkdir -p "$folder_genomics/assembly/$genome_ids[$i]"
    zsh 01a-filter_reads.sh \
        $raw_reads \
        $filtered_reads

    # Extract the filtered reads to a txt file for plotting
    zsh 01b-extract_reads.sh \
        $filtered_reads \
        "$folder_genomics/assembly/$genome_ids[$i]/filtered_reads.txt" 
    
    # Plot the read data
    Rscript 01c-plot_reads.R \
        "$folder_genomics/assembly/$genome_ids[$i]/filtered_reads.txt"  \
        "$folder_genomics/assembly/$genome_ids[$i]/filtered_reads_qc.png" 
done



