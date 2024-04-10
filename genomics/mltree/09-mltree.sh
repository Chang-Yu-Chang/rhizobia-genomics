#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script contrusted ML tree

cd $folder_shell
mkdir -p $folder_genomics/mltree


mamba activate iqtree

# Test on one gene
aln='aarA'
mkdir -p $folder_genomics/mltree/$aln
iqtree \
    -nt AUTO \
    -s $folder_genomics/pangenome/panaroo/aligned_gene_sequences/$aln.aln.fas \
    --prefix $folder_genomics/mltree/$aln/aln

# All core genes
mkdir -p $folder_genomics/mltree/core_b
iqtree \
    -nt AUTO \
    -B 1000 \
    -s $folder_genomics/pangenome/panaroo/core_gene_alignment.aln \
    --prefix $folder_genomics/mltree/core_b/aln
