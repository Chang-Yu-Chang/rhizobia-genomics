#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# Single copy core genes
# Run IQ-TREE on the combined alignment file
mamba activate iqtree

mkdir -p $folder_genomics/pangenome/trees/combined_sccg
iqtree \
    -nt AUTO \
    -B 1000 \
    -s $folder_genomics/pangenome/trees/combined_sccg_alignment.fas \
    --prefix $folder_genomics/pangenome/trees/combined_sccg/combined_sccg
