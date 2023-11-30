#!/usr/bin/env zshs
source ~/.zshrc

# This downsamples the reads to 250 Mbp
# $1: filtered_reads in fastq
# $2: downsampled_reads in fastq


conda activate
mamba activate filtlong

filtlong --target_bases 250000000 $1 | gzip > $2
# `--target_bases 250000000` remove the worst reads unitl only 250 Mbp remain
