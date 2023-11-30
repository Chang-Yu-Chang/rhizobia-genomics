#!/usr/bin/env zsh
source ~/.zshrc

# This indexes the references genomes
# $1: reference genome in mmi
# $2: reference genome in fasta

conda activate
mamba activate minimap2

minimap2 -d $1 $2

