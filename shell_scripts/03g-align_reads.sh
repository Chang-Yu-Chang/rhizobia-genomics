#!/usr/bin/env zsh
source ~/.zshrc

# This aligns the filtered reads to the reference genome
# $1: reference genome in fasta
# $2: target reads in fastq
# $3: target genome sam file

conda activate
mamba activate minimap2

minimap2 -ax map-ont $1 $2 > $3
