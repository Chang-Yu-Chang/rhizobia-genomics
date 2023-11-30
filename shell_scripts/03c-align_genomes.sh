#!/usr/bin/env zsh
source ~/.zshrc

# This aligns the genome to the reference
# $1: reference genome in mmi
# $2: target genome
# $3: target genome sam file

conda activate
mamba activate minimap2

minimap2 -ax map-ont $1 "$folder_anvio/01-fasta/$i.fasta" > "$folder_anvio/06-alignment/sam/$i.sam"
