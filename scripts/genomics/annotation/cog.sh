#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# Classify prokaryote protein sequences into COG functional category.
mkdir -p $folder_genomics/cog
mamba activate cogclassifier

for i in {1..38}
do
    echo $genome_ids[$i]
    COGclassifier \
        -i $folder_genomics/faa/$genome_ids[$i].faa \
        -o $folder_genomics/cog/$genome_ids[$i]/
done
