#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script constructs ML trees using iqtree

cd $folder_shell
mkdir -p $folder_genomics/mltree

mamba activate iqtree

mkdir -p $folder_genomics/mltree/core_b
iqtree \
    -nt AUTO \
    -B 1000 \
    -s $folder_genomics/pangenome/core_gene_alignment.aln \
    --prefix $folder_genomics/mltree/core_b/aln

# -nt AUTO will tell IQ-TREE to automatically determine the best number of cores given the current data and computer.
# -B Replicates for ultrafast bootstrap (>=1000)
