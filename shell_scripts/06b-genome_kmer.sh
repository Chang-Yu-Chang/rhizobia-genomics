#!/usr/bin/env zsh
source ~/.zshrc

# This creates sourmash kmer signature from genome
# $1: genome in fasta
# $2: output kmer signature

mamba activate sourmash

sourmash sketch dna -p k=31,scaled=1000 -o $2 $1
