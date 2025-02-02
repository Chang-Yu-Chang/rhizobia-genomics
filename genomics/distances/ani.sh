#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

cd $folder_shell

# This script computes pairwise genomic ANI

mkdir -p $folder_genomics/distances/ani

## Create a list of genome fasta files
for i in {1..38}
do
    echo -e $folder_genomics/fasta/genomes/$genome_ids[$i].fasta
done >| $folder_genomics/distances/ani/list_genomes.txt

# echo $folder_genomics/fasta/genomes/em1021.fasta >> $folder_genomics/ani/list_genomes.txt
# echo $folder_genomics/fasta/genomes/em1022.fasta >> $folder_genomics/ani/list_genomes.txt
# echo $folder_genomics/fasta/genomes/usda1106.fasta >> $folder_genomics/ani/list_genomes.txt
# echo $folder_genomics/fasta/genomes/wsm419.fasta >> $folder_genomics/ani/list_genomes.txt
# echo $folder_genomics/fasta/genomes/casidaa.fasta >> $folder_genomics/ani/list_genomes.txt

## Compute ani
zsh 06a-fastani.sh \
    $folder_genomics/distances/ani/list_genomes.txt \
    $folder_genomics/distances/ani/ani_genomes.txt
