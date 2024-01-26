#!/usr/bin/env zsh
source ~/.zshrc

# This extracts the rRNA sequences from the genome fasta and blast
# $1: $folder_genomics

# GCF_002197065.1 is E. meliloti usda1106
# GCF_000006965.1 is E. meliloti em1021
# GCF_013315775.1 is E. meliloti em1022
# GCF_000017145.1 is E. medicae wsm419

mamba activate blast

mkdir -p "$1/blast_db"
cat \
    "$1/genomes/em1021.fasta" \
    "$1/genomes/em1022.fasta" \
    "$1/genomes/usda1106.fasta" \
    "$1/genomes/wsm419.fasta" \
    > "$1/blast_db/genomes.fasta"

makeblastdb -in "$1/blast_db/genomes.fasta" -dbtype nucl -out "$1/blast_db/genomes"











