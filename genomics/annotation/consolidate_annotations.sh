#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script consolidates annotation files into one folder depending on the file type

## gff
mkdir -p $folder_genomics/gff/genomes

for i in {1..32}; do
    cp $folder_genomics/annotation/genomes/$genome_ids[$i]/annotated.gff $folder_genomics/gff/genomes/$genome_ids[$i].gff
done

for ref in em1021 em1022 usda1106 wsm419 casidaa; do
    cp $folder_genomics/annotation/genomes/$ref/annotated.gff $folder_genomics/gff/genomes/$ref.gff
done

## faa
mkdir -p $folder_genomics/faa/genomes

for i in {1..32}; do
    cp $folder_genomics/annotation/genomes/$genome_ids[$i]/annotated.faa $folder_genomics/faa/genomes/$genome_ids[$i].faa
done

for ref in em1021 em1022 usda1106 wsm419 casidaa; do
    cp $folder_genomics/annotation/genomes/$ref/annotated.faa $folder_genomics/faa/genomes/$ref.faa
done
