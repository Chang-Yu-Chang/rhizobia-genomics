#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

mamba activate biopython

set_name='elev_med'
mkdir -p $folder_data/phylogenomics_analysis/trees/$set_name
echo -e "g4\ng5\ng6\ng8\ng9\ng11\ng13\ng16\ng17\ng19" > $folder_data/phylogenomics_analysis/trees/elev_med/list_genomes.csv
# single copy gene list
cp $folder_data/genomics_analysis/gene_content/$set_name/list_sccg.csv $folder_data/phylogenomics_analysis/trees/$set_name/list_sccg.csv

# Concatenate
python concatenate_alignment.py \
    $folder_genomics/pangenome/$set_name/aligned_gene_sequences \
    $folder_genomics/pangenome/$set_name/combined_sccg_alignment.fas \
    $folder_data/phylogenomics_analysis/trees/$set_name/list_sccg.csv \
    $folder_data/phylogenomics_analysis/trees/$set_name/list_genomes.csv
cp $folder_genomics/pangenome/$set_name/combined_sccg_alignment.fas $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg_alignment.fas


set_name='urbn_mel'
mkdir -p $folder_data/phylogenomics_analysis/trees/urbn_mel
echo -e "g21\ng22\ng23\ng24\ng25\ng26\ng27\ng31\ng32\ng33\ng34\ng35\ng36\ng37\ng39\ng41\ng42\ng43\ng44\ng45" > $folder_data/phylogenomics_analysis/trees/urbn_mel/list_genomes.csv
# single copy gene list
cp $folder_data/genomics_analysis/gene_content/$set_name/list_sccg.csv $folder_data/phylogenomics_analysis/trees/$set_name/list_sccg.csv

# Concatenate
python concatenate_alignment.py \
    $folder_genomics/pangenome/$set_name/aligned_gene_sequences \
    $folder_genomics/pangenome/$set_name/combined_sccg_alignment.fas \
    $folder_data/phylogenomics_analysis/trees/$set_name/list_sccg.csv \
    $folder_data/phylogenomics_analysis/trees/$set_name/list_genomes.csv
cp $folder_genomics/pangenome/$set_name/combined_sccg_alignment.fas $folder_data/phylogenomics_analysis/trees/$set_name/combined_sccg_alignment.fas
