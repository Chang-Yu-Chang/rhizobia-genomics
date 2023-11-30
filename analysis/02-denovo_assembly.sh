#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script performs de novo genome assembly

cd $folder_shell
echo "02-denovo_assembly"

for i in {1..19}
do
    echo "$folder_raw/$batch_names[$i]/$sample_ids[$i]"
    raw_reads="$folder_raw/$batch_names[$i]/$sample_ids[$i]/reads/raw_reads.fastq.gz"
    filtered_reads="$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.fastq.gz"
    downsampled_reads="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/downsampled_reads.fastq"
    assembly_prefix="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/miniasm/draft_assembly"
    draft_genome="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/miniasm/draft_genome.fasta"
    downsampled_reads2="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/downsampled_reads2.fastq"

    # Downsample reads
    zsh 02a-downsample_reads.sh \
        $filtered_reads \
        $downsampled_reads

    # Draft genome via miniasm
    mkdir -p "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/miniasm"
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
    mkdir -p "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/flye"
    zsh 02d-flye_assembly.sh \
        $downsampled_reads2 \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/flye"

    # Polish the assembled genome via medaka
    mkdir -p "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/medaka"
    zsh 02e-medaka_clean.sh \
        $raw_reads \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/flye/assembly.fasta" \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/medaka"

    # Clean up consensus name
    zsh 02f-clean_names.sh \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/medaka/consensus.fasta" \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/cleaned_names.txt"

    # Check assembly quality via quast
    mkdir -p "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/quast"
    zsh 02g-quast.sh \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"\
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/quast"

    # Check assembly quality via busco
    mkdir -p "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/busco"
    zsh 02h-busco.sh \
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"\
        "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/busco"


done





