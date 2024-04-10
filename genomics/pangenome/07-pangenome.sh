#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
mkdir -p $folder_genomics/pangenome

# 1. Whole genomes
mkdir -p $folder_genomics/pangenome/genomes

# Create a list of gff
for i in {1..22} {24..32}; do; echo -e $folder_genomics/gff/$genome_ids[$i].gff
done >| $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/em1021.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/em1022.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/usda1106.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/wsm419.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/casidaa.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt

zsh 07b-panaroo.sh \
    $folder_genomics/pangenome/genomes \
    $folder_genomics/pangenome/genomes/list_gffs.txt

# 2. contigs
mkdir -p $folder_genomics/pangenome/contigs

# Create a list of gff
for i in {1..22} {24..32}; do; echo -e $folder_genomics/gff/$genome_ids[$i].gff
done >| $folder_genomics/pangenome/contigs/list_gffs.txt
echo $folder_genomics/gff/em1021.gff >> $folder_genomics/pangenome/contigs/list_gffs.txt
echo $folder_genomics/gff/em1022.gff >> $folder_genomics/pangenome/contigs/list_gffs.txt
echo $folder_genomics/gff/usda1106.gff >> $folder_genomics/pangenome/contigs/list_gffs.txt
echo $folder_genomics/gff/wsm419.gff >> $folder_genomics/pangenome/contigs/list_gffs.txt
echo $folder_genomics/gff/casidaa.gff >> $folder_genomics/pangenome/contigs/list_gffs.txt

zsh 07b-panaroo.sh \
    $folder_genomics/pangenome/contigs \
    $folder_genomics/pangenome/contigs/list_gffs.txt




# 1. Roary
# mkdir -p $folder_genomics/pangenome/roary
# zsh 07a-roary.sh \
#     $folder_genomics/pangenome/roary \
#     $folder_genomics/gff



# mamba activate panaroo
# panaroo-integrate -d $folder_genomics/pangenome/panaroo/ -i $folder_genomics/gff/casidaa.gff -t 24 -o $folder_genomics/pangenome/updated

# Infinite Many Genes model
#mamba activate panaroo
#panaroo-fmg --tree dated_phylogeny.newick --pa gene_presence_absence_renamed.Rtab -o fmg_results.txt
# # Subset the core gene sequence
# cd $folder_genomics/pangenome/panaroo/aligned_gene_sequences
# ls *.aln.fas | sed 's/\.aln\.fas$//' > ../list_core_gene.txt
#
# mamba activate seqtk
# seqtk subseq pan_genome_reference.fa list_core_gene.txt > core_gene.fasta
#
# # Call variants from aligned genes
# #zsh 07c-core_genes.sh
# cd L
