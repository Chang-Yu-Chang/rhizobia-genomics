#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# Single copy core genes
# Run IQ-TREE on the combined alignment file
mamba activate iqtree

mkdir -p $folder_data/phylogenomics_analysis/trees/combined_sccg
iqtree \
    -nt AUTO \
    -B 1000 \
    -s $folder_data/phylogenomics_analysis/trees/combined_sccg_alignment.fas \
    --prefix $folder_data/phylogenomics_analysis/trees/combined_sccg/combined_sccg
