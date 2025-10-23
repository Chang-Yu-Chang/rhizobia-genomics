#!/usr/bin/env zsh
# 1. `concatenate_alignment.sh` concatenates the single-copy core-gene alignment into one alignment file. It also copies the aligment file to the folder tree/
# 1a. `concatenate_alignment.py` is the command line tool for 1.

source ~/.zshrc
source ../env_vars.sh

mamba activate biopython
mkdir -p $folder_genomics/pangenome/trees/
echo $genome_ids | tr ' ' '\n' > $folder_genomics/pangenome/trees/list_genomes.csv

# Concatenate
python concatenate_alignment.py \
    $folder_genomics/pangenome/gene_content/aligned_gene_sequences \
    $folder_genomics/pangenome/gene_content/combined_sccg_alignment.fas \
    $folder_genomics/pangenome/gene_content/list_sccg.csv \
    $folder_genomics/pangenome/gene_content/list_genomes.csv

# Copy this alignment for to the tree folders
cp $folder_genomics/pangenome/panaroo/combined_sccg_alignment.fas $folder_genomics/trees/combined_sccg_alignment.fas
