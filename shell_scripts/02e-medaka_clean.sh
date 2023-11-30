#!/usr/bin/env zsh
source ~/.zshrc

# This polishes flye genome via medaka
# $1: raw_reads=$1 in fastq.gz
# $2: draft_assembly in fasta
# $3: medaka_folder

conda activate
mamba activate medaka

medaka_consensus -i $1 -d $2 -o $3 -t 10 -m r941_min_high_g303
# `-i "$raw_reads"` fastx input basecalls (required).
# `-d "$draft_assembly"` fasta input assembly (required).
# `-o "$medaka_folder"` output folder (default: medaka).
# `-t 10` number of threads with which to create features (default: 1).
# `-m r941_min_high_g303` medaka model, (default: r941_min_high_g360).
