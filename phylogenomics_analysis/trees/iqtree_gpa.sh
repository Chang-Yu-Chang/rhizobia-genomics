#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script constructs ML trees using iqtree
mamba activate iqtree

for phy in gpa_genomes gpa_chrom gpa_psyma gpa_psymb; do;
    iqtree \
        -nt AUTO \
        -bb 1000 -nt AUTO -m MK+R+FO -redo -st MORPH \
        -s $folder_data/phylogenomics_analysis/trees/mltree/$phy/$phy.phy \
        --prefix $folder_data/phylogenomics_analysis/trees/mltree/$phy/$phy
done
