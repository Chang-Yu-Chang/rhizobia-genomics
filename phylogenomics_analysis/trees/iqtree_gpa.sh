#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script constructs ML trees using iqtree

cd $folder_shell
mkdir -p $folder_genomics/mltree

# Create the isolate single-copy core gene tree
mkdir -p $folder_genomics/mltree/isolates_gpa

# Run IQ-TREE on the combined alignment file
mamba activate iqtree

iqtree \
    -nt AUTO \
    -bb 1000 -nt AUTO -m MK+R+FO -redo -st MORPH \
    -s $folder_genomics/mltree/isolates_gpa/gpa_binary.phy \
    --prefix $folder_data/phylogenomics_analysis/trees/mltree/isolates_gpa/aln

# -nt AUTO will tell IQ-TREE to automatically determine the best number of cores given the current data and computer.
# -B Replicates for ultrafast bootstrap (>=1000)
