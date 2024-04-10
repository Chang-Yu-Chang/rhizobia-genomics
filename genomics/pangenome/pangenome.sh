#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
mkdir -p $folder_genomics/pangenome

mkdir -p $folder_genomics/pangenome/genomes

# Create a list of gff. Skip g28
for i in {1..22} {24..32}
do
    echo -e $folder_genomics/gff/$genome_ids[$i].gff
done >| $folder_genomics/pangenome/list_gffs.txt

echo $folder_genomics/gff/em1021.gff >> $folder_genomics/pangenome/list_gffs.txt
echo $folder_genomics/gff/em1022.gff >> $folder_genomics/pangenome/list_gffs.txt
echo $folder_genomics/gff/usda1106.gff >> $folder_genomics/pangenome/list_gffs.txt
echo $folder_genomics/gff/wsm419.gff >> $folder_genomics/pangenome/list_gffs.txt
echo $folder_genomics/gff/casidaa.gff >> $folder_genomics/pangenome/list_gffs.txt

zsh 07b-panaroo.sh \
    $folder_genomics/pangenome/genomes \
    $folder_genomics/pangenome/list_gffs.txt
