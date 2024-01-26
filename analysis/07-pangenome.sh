#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
mkdir -p $folder_genomics/pangenome

# 1. Roary
mkdir -p $folder_genomics/pangenome/roary
zsh 07a-roary.sh \
    $folder_genomics/pangenome/roary \
    $folder_genomics/gff

# 2. Panaroo
mkdir -p $folder_genomics/pangenome/panaroo

# Create a list of gff
for i in {1..38}; do; echo -e $folder_genomics/gff/$genome_ids[$i].gff
done >| $folder_genomics/pangenome/panaroo/list_gffs.txt

# for i in {1..5}; do; echo -e $folder_genomics/gff/$genome_ids[$i].gff
# done >| $folder_genomics/pangenome/panaroo/list_gffs.txt

# Generate gene presence-absence table
# It also performs multiple sequence alignment of core genes using MAFFT
zsh 07b-panaroo.sh \
    $folder_genomics/pangenome/panaroo \
    $folder_genomics/pangenome/panaroo/list_gffs.txt

# Infinite Many Genes model
#mamba activate panaroo
#panaroo-fmg --tree dated_phylogeny.newick --pa gene_presence_absence_renamed.Rtab -o fmg_results.txt

# Subset the core gene sequence
cd $folder_genomics/pangenome/panaroo/aligned_gene_sequences
ls *.aln.fas | sed 's/\.aln\.fas$//' > ../list_core_gene.txt

mamba activate seqtk
seqtk subseq pan_genome_reference.fa list_core_gene.txt > core_gene.fasta

# Call variants from aligned genes
zsh 07c-core_genes.sh

