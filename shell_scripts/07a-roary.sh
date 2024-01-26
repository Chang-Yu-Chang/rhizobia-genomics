#!/usr/bin/env zsh
source ~/.zshrc

# This uses roary to analyze pangenome
# $1: output roary folder
# $2: prokka annotation gffs

conda activate
mamba activate roary

roary -f $1 -e -n -i 95 -cd 99 -s $2
# `-f STR` output directory
# `-e` create a multiFASTA alignment of core genes using PRANK
# `-n` fast core gene alignment with MAFFT
# `-i 95` minimum percentage identity for blastq [95]
# `-cd FLOAT` percentage of isolates in a gene must be in to be core [99]
# `-r` create R plots
# `-s` dont split paralogs
# `*.gff` input
