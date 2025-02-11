#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assembles the genomes from raw reads following the plasmidsaurus pipeline

cd $folder_shell

for i in {1..38}
do
    echo $genome_ids[$i]
    dir=$folder_genomics/assembly/$genome_ids[$i]
    raw_reads=$folder_genomics/raw_reads/$genome_ids[$i].fastq.gz
    filtered_reads=$dir/filtered_reads.fastq.gz
    downsampled_reads=$dir/downsampled_reads.fastq
    assembly_prefix=$dir/miniasm/draft_assembly
    draft_genome=$dir/miniasm/draft_genome.fasta
    downsampled_reads2=$dir/downsampled_reads2.fastq

    # Downsample reads
    zsh 02a-downsample_reads.sh \
        $filtered_reads \
        $downsampled_reads

    # Draft genome via miniasm
    mkdir -p $dir/miniasm
    zsh 02b-draft_genome.sh \
        $downsampled_reads \
        $assembly_prefix \
        $draft_genome

    # Downsample reads again
    genome_size=$(grep -v ">" $draft_genome | grep -E -o "G|C|T|A|N" | wc -l)
    zsh 02c-downsample_reads2.sh \
        $genome_size \
        $draft_genome \
        $filtered_reads \
        $downsampled_reads2

    # Assemble genome via flye
    mkdir -p $dir/flye
    zsh 02d-flye_assembly.sh \
        $downsampled_reads2 \
        $dir/flye

    # Polish the assembled genome via medaka
    mkdir -p $dir/medaka
    zsh 02e-medaka_clean.sh \
        $raw_reads \
        $dir/flye/assembly.fasta \
        $dir/medaka

done





