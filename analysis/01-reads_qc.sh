#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script performs quality control on the raw nanapore whole-genome long reads

cd $folder_shell

for i in {1..19}
do
    echo "$folder_raw/$batch_ids[$i]/$sample_ids[$i]"
    raw_reads="$folder_raw/$batch_ids[$i]/$sample_ids[$i]/reads/raw_reads.fastq.gz"
    filtered_reads="$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.fastq.gz"

    # Filter the worst 5% reads via filtlong
    zsh 01-filter_reads.sh \
        $raw_reads \
        $filtered_reads

    # Extract the filtered reads to a txt file
    zsh 01a-extract_reads.sh \
        $filtered_reads \
        "$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.txt"

    # Plot the read data
    Rscript 01b-plot_reads.R \
        "$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.txt" \
        "$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads_qc.png"
done



