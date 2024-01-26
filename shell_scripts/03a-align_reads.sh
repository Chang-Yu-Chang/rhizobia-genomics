#!/usr/bin/env zsh
source ~/.zshrc

# This aligns the filtered reads to the reference genome and calls SNPs
# $1: reference genome in fasta
# $2: target reads in fastq
# $3: snippy output folder

conda activate
mamba activate snippy

snippy --ref $1 --ctgs $2 --outdir $3 --force


