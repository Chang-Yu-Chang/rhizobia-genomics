#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh


mamba activate biopython

set_name='elev_med'
echo -e "g4\ng5\ng6\ng8\ng9\ng11\ng13\ng16\ng17\ng19" > $folder_data/phylogenomics_analysis/replicon_trees/elev_med/list_genomes.csv

for replicon in 'chromosome' 'pSymA' 'pSymB' "pAcce"; do
    mkdir -p $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon

    # Concatenate
    python concatenate_alignment.py \
        $folder_genomics/pangenome/$set_name/aligned_gene_sequences \
        $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg_alignment.fas \
        $folder_data/phylogenomics_analysis/replicon_trees/$set_name/list_sccg_$replicon.csv \
        $folder_data/phylogenomics_analysis/replicon_trees/$set_name/list_genomes.csv
done


set_name='urbn_mel'
echo -e "g21\ng22\ng23\ng24\ng25\ng26\ng27\ng31\ng32\ng33\ng34\ng35\ng36\ng37\ng39\ng41\ng43" > $folder_data/phylogenomics_analysis/replicon_trees/urbn_mel/list_genomes.csv

for replicon in 'chromosome' 'pSymA' 'pSymB'; do #  no single copy core on pAcce
    mkdir -p $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon

    # Concatenate
    python concatenate_alignment.py \
        $folder_genomics/pangenome/$set_name/aligned_gene_sequences \
        $folder_data/phylogenomics_analysis/replicon_trees/$set_name/$replicon/combined_sccg_alignment.fas \
        $folder_data/phylogenomics_analysis/replicon_trees/$set_name/list_sccg_$replicon.csv \
        $folder_data/phylogenomics_analysis/replicon_trees/$set_name/list_genomes.csv
done
