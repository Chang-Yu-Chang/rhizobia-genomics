#!/usr/bin/env zsh
source ~/.zshrc

# This selects for high quality ONT reads and assemble again
# $1: downsampled_reads2 in fastq
# $2: flye folder

mamba activate flye

flye --nano-raw $1 --out-dir $2 --threads 20 --genome-size '7m'
# `--meta` metagenome / uneven coverage mode
# --nano-raw path [path ...] ONT regular reads, pre-Guppy5 (<20% error)
# `--nano-corr "$downsampled_reads2"` ONT reads that were corrected with other methods (<3% error)
# `--out-dir "$flye_folder"` Output directory
# `--threads 20`  number of parallel threads
# `--genome-size 7m` estimated genome size 7 Mb
# The --min-overlap parameter controls the minimum overlap length between reads for them to be merged. Increasing this value can sometimes help the assembler be more conservative about splitting contigs. However, this can also reduce sensitivity, so you need to balance it carefully.`
