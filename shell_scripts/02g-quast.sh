#!/usr/bin/env zsh
source ~/.zshrc

# This checks genome quality
# $1: genome in fasta
# $2: quast folder

conda activate
mamba activate quast

quast $1 -o $2
# `-o` output directory
