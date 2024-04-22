#!/usr/bin/env zsh
source ~/.zshrc

# This calculates the genomic ANI between a list of genomes
# $1: list of genomes for ani in txt
# $2: fasani output in out

mamba activate fastani

fastANI -t 10 \
    --ql $1 \
    --rl $1 \
    -o $2
# `--ql` list of names of query sequences in fasta
# `--rl` list of names of reference sequences in fasta
