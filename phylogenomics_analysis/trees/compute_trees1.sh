#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# Single copy core genes
# Run IQ-TREE on the combined alignment file
mamba activate iqtree

set_name='elev_med'
mkdir -p $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg
iqtree \
-nt AUTO \
-B 1000 \
-s $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg_alignment.fas \
--prefix $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg/combined_sccg

set_name='urbn_mel'
mkdir -p $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg
iqtree \
-nt AUTO \
-B 1000 \
-s $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg_alignment.fas \
--prefix $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg/combined_sccg
