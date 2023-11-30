#!/usr/bin/env zsh
source ~/.zshrc

# This allows genome comparison
# $1: genome in fasta
# $2: sourmash_folder
# $3: gtdb_db

conda activate
mamba activate sourmash

sourmash_sig="$2/consensus.sig"
gathered_csv="$2/gathered.csv"

# Create a signature for this genome
sourmash sketch dna --check-sequence -f -p scaled=1000,k=31 $1 -o "$2/consensus.sig"
# `sketch dna` command reads in DNA sequences and outpouts DNA sketches
# `-p scaled=1000,k=31` signature parameters to use
# `-o` output computed signatures to this directory

# We recommend using the Zipfile databases for sourmash gather and the SBT databases for sourmash search.
sourmash gather $2/consensus.sig $3 --save-matches $2/matches.zip



