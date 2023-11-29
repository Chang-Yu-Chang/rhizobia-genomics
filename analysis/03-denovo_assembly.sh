#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script performs de novo genome assembly

cd $folder_shell
echo "03-denovo_assembly"

# 2. Draft genome via miniasm
echo "2-miniasm"

downsample_reads="$folder_genomes/$sample_id/02-miniasm/02-downsampled_reads.fastq"
assemble_paf="$folder_genomes/$sample_id/02-miniasm/02-miniasm_pilot.paf.gz"
assemble_gfa="$folder_genomes/$sample_id/02-miniasm/02-miniasm_pilot.gfa"
assemble_fq="$folder_genomes/$sample_id//02-miniasm/02-miniasm_pilot.fasta"

zsh 32-miniasm.sh $filtered_reads $downsample_reads $assemble_paf $assemble_gfa $assemble_fq &> $log02

# 3. Assemble genome via flye
echo "3-flye"
filtered_reads="$folder_genomes/$sample_id/01-filtlong/01-filtered_reads.fastq.gz"
assemble_fq="$folder_genomes/$sample_id/02-miniasm/02-miniasm_pilot.fasta"
downsampled_reads2="$folder_genomes/$sample_id/03-flye/03-downsampled2_reads.fastq"
genome_size=$(grep -v ">" "$assemble_fq" | grep -E -o "G|C|T|A|N" | wc -l)
flye_folder="$folder_genomes/$sample_id/03-flye/"

zsh 33-flye.sh $filtered_reads $assemble_fq $downsampled_reads2 $genome_size $flye_folder &> $log03

# 4. Polish genome via medaka
echo "4-medaka"
raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
draft_assembly="$folder_genomes/$sample_id/03-flye/assembly.fasta"
medaka_folder="$folder_genomes/$sample_id/04-medaka"
