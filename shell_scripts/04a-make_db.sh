#!/usr/bin/env zsh
source ~/.zshrc

# This extracts the rRNA sequences from the genome fasta and blast
# $1: $folder_genomics
# $2: list of genomes

mkdir -p $1/blast_db
list_strains=("${(@f)$(cat $2 | cut -d ',' -f 3)}")
n_strains="${#list_strains}"
for i in {1..$n_strains}; do
    cat $1/genomes/$list_strains[$i].fasta
done |> $1/blast_db/genomes.fasta

mamba activate blast

makeblastdb -in $1/blast_db/genomes.fasta -dbtype nucl -out $1/blast_db/genomes











