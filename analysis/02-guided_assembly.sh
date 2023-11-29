#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script performs reference-guided genome assembly

cd $folder_shell
echo "02-guided_assembly"

zsh 02a-align_genomes.sh  \
    "$folder_genomics/reference/usda1106.mmi" \
    $sample_id


