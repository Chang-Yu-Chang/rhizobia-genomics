#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# Create the isolate single-copy core gene tree
mkdir -p $folder_data/phylogenomics_analysis/trees/mltree/seq_core

# Run IQ-TREE on the combined alignment file
mamba activate iqtree

iqtree \
    -nt AUTO \
    -B 1000 \
    -s $folder_genomics/pangenome/isolates/core_gene_alignment.aln \
    --prefix $folder_data/phylogenomics_analysis/trees/mltree/seq_core/seq_core

# -nt AUTO will tell IQ-TREE to automatically determine the best number of cores given the current data and computer.
# -B Replicates for ultrafast bootstrap (>=1000)
