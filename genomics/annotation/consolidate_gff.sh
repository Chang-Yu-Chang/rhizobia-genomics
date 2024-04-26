#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script consolidates genome gff into one folder

mkdir -p $folder_genomics/gff/genomes

for i in {1..32}; do
    cp $folder_genomics/annotation/genomes/$genome_ids[$i]/annotated.gff $folder_genomics/gff/genomes/$genome_ids[$i].gff
done

for ref in em1021 em1022 usda1106 wsm419 casidaa; do
    cp $folder_genomics/annotation/genomes/$ref/annotated.gff $folder_genomics/gff/genomes/$ref.gff
done
