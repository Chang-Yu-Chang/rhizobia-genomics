#!/usr/bin/env zsh
source ~/.zshrc

# This selects for high quality ONT reads and assemble again
# $1: downsampled_reads2 in fastq
# $2: flye folder

mamba activate flye

flye --meta --nano-corr $1 --out-dir $2 --threads 20
# `--meta` metagenome / uneven coverage mode
# `--nano-corr "$downsampled_reads2"` ONT reads that were corrected with other methods (<3% error)
# `--out-dir "$flye_folder"` Output directory
# `--threads 20`  number of parallel threads


