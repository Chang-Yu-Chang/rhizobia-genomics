#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script finds the orthologous genes
mamba activate orthofinder

orthofinder \
    -t 10 \
    -f $folder_data/genomics/faa/genomes \
    -o $folder_data/phylogenomics_analysis/orthofinder/genomes
