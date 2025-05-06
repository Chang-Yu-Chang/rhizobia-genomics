#!/usr/bin/env zsh
source ~/.zshrc

# This compares the kmer signatures
# $1: list of signatures
# $2: output distance matrix

mamba activate sourmash

# Create sourmash kmer signatures
sourmash compare \
    --from-file $1 \
    --ksize 31 \
    --distance-matrix \
    --csv $2
