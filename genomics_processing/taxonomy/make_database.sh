#!/usr/bin/env zsh
source ~/.zshrc
source ../00-env_vars.sh

# This script makes a custom BLAST genome database

mkdir -p $folder_genomics/blast_db

list_strains=("${(@f)$(cat $folder_data/raw/ensifer_ncbi.csv | cut -d ',' -f 3)}")
n_strains="${#list_strains}"

for i in {1..$n_strains}; do
    cat $folder_genomics/fasta/genomes/$list_strains[$i].fasta
done |> $folder_genomics/blast_db/genomes.fasta

mamba activate blast

makeblastdb -in $folder_genomics/blast_db/genomes.fasta -dbtype nucl -out $folder_genomics/blast_db/genomes

