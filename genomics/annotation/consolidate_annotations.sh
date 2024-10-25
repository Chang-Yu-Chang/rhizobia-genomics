#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script consolidates annotation files into one folder depending on the file type

## gff
mkdir -p $folder_genomics/gff

for i in {1..38}; do
    cp $folder_genomics/annotation/$genome_ids[$i]/annotated.gff $folder_genomics/gff/$genome_ids[$i].gff
done

# for ref in em1021 em1022 usda1106 wsm419 casidaa; do
for ref in em1021 wsm419 ngr234; do
    cp $folder_genomics/annotation/$ref/annotated.gff $folder_genomics/gff/$ref.gff
done

## faa
mkdir -p $folder_genomics/faa

for i in {1..38}; do
    cp $folder_genomics/annotation/$genome_ids[$i]/annotated.faa $folder_genomics/faa/$genome_ids[$i].faa
done

#for ref in em1021 em1022 usda1106 wsm419 casidaa; do
for ref in em1021 wsm419 ngr234; do
    cp $folder_genomics/annotation/$ref/annotated.faa $folder_genomics/faa/$ref.faa
done
