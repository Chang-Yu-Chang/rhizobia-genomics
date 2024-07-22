#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script constructs ML trees using iqtree

cd $folder_shell
mkdir -p $folder_genomics/mltree

# Create the isolate single-copy core gene tree
mkdir -p $folder_genomics/mltree/isolates_core_b

# Run IQ-TREE on the combined alignment file
mamba activate iqtree

iqtree \
    -nt AUTO \
    -B 1000 \
    -s $folder_genomics/pangenome/isolates/core_gene_alignment.aln \
    --prefix $folder_data/phylogenomics_analysis/treesmltree/isolates_core_b/aln

# -nt AUTO will tell IQ-TREE to automatically determine the best number of cores given the current data and computer.
# -B Replicates for ultrafast bootstrap (>=1000)
