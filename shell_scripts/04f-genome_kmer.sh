#!/usr/bin/env zsh
source ~/.zshrc

# This blast the whole genome across the database
# $1: genome in fasta
# $2: output kmer signature

conda activate
mamba activate sourmash

# Create sourmash kmer signatyre
sourmash compute -k 31 --scaled 1000 -o $2 $1
