#!/usr/bin/env zsh
source ~/.zshrc
source ../../env_vars.sh

# This script finds the orthologous genes
mamba activate orthofinder

orthofinder \
    -f $folder_genomics/faa/genomes \
    -o $folder_data/phylogenomics_analysis/orthofinder/genomes
