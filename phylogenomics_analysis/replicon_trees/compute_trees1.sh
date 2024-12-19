#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# Single copy core genes
# Run IQ-TREE on the combined alignment file
mamba activate iqtree

set_name='elev_med'

for replicon in 'chromosome' 'pSymA' 'pSymB' "pAcce"; do
    mkdir -p $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg

    iqtree \
        -nt AUTO \
        -B 1000 \
        -s $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg_alignment.fas \
        --prefix $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg/combined_sccg
done


set_name='urbn_mel'
for replicon in 'chromosome' 'pSymA' 'pSymB'; do
    mkdir -p $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg

    iqtree \
        -nt AUTO \
        -B 1000 \
        -s $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg_alignment.fas \
        --prefix $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg/combined_sccg
done
