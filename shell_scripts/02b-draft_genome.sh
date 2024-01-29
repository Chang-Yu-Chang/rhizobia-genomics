#!/usr/bin/env zsh
source ~/.zshrc

# This drafts a genome
# $1: downsampled_reads in fastq
# $2: assembly_prefix
# $3: draft_genome in fasta

mamba activate miniasm

# Overlap raw reads using minimap
minimap2 -x ava-ont -t8 $1 $1 | gzip -1 > "$2.paf.gz"
# `-x ava-ont` present Nanopore read overlap
# `-t8` use 8 threads
# `gzip -1` fastest compression

# Assemble reads using miniams
miniasm -f $1 "$2.paf.gz"> "$2.gfa"
# `-f` read sequences

# Convert assembly graph to fastq
any2fasta "$2.gfa" > $3
