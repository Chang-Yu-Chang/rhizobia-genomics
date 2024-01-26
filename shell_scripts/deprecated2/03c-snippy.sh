#!/usr/bin/env zsh
source ~/.zshrc

# This calls SNPs via snippy
# $1: reference genome in gbff
# $2: raw reads in fasta
# $3: snippy output

conda activate
mamba activate snippy

# For snippy, I dont need to index my fasta. I can use the fasta file directly
snippy --ref $1 --ctgs $2 --outdir $3 --force


