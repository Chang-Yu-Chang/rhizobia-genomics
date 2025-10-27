#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script annotates genomes

cd $folder_shell

mkdir -p $folder_genomics/annotation

## Assembled genomes
for i in {1..38}
do
    echo $genome_ids[$i]
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$genome_ids[$i]/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/$genome_ids[$i].fasta \
        $folder_genomics/annotation/$genome_ids[$i]/
done

## gff
mkdir -p $folder_genomics/gff

for i in {1..38}; do
    cp $folder_genomics/annotation/$genome_ids[$i]/annotated.gff $folder_genomics/gff/$genome_ids[$i].gff
done
