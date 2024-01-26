#!/usr/bin/env zsh
source ~/.zshrc

# This cleans the assembled genome fasta to have simplified deflines
# $1: medaka_consensus in fasta
# $2: genome name in fasta
# $3 cleaned names mapping file in txt

conda activate
mamba activate anvio-8


anvi-script-reformat-fasta \
    $1 -o $2 \
    --simplify-names \
    --report-file $3 \
    --min-len 500
# --simplify-names edit deflines to make sure they have simple names. This is important for BAM
# --report-file reports the changes to deflines
# --min-len 500 minimum length of contigs to keep



