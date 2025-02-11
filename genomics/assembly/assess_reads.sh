#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assesses the quality of the raw nanapore whole-genome long reads

cd $folder_shell
mkdir -p $folder_genomics/assembly

for i in {1..38}
do
    echo $genome_ids[$i]
    raw_reads=$folder_genomics/raw_reads/$genome_ids[$i].fastq.gz
    filtered_reads=$folder_genomics/assembly/$genome_ids[$i]/filtered_reads.fastq.gz
    mkdir -p $folder_genomics/assembly/$genome_ids[$i]

    # Remove the worst 5% reads via filtlong
    zsh 01a-filter_reads.sh \
        $raw_reads \
        $filtered_reads

    # # Extract the filtered reads to a txt file for plotting
    zsh 01b-extract_reads.sh \
        $filtered_reads \
        $folder_genomics/assembly/$genome_ids[$i]/filtered_reads.txt

    # # Plot the filtered reads
    # Rscript 01c-plot_reads.R \
    #     $folder_genomics/assembly/$genome_ids[$i]/filtered_reads.txt \
    #     $folder_genomics/assembly/$genome_ids[$i]/filtered_reads_qc.png
done
