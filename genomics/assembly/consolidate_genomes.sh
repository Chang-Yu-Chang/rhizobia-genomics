#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script consolidates genome fasta files into one folderx

for i in {1..38}
do
    dir=$folder_genomics/assembly/$genome_ids[$i]
    cp $dir/medaka/consensus.fasta $folder_genomics/fasta/genomes/$genome_ids[$i].fasta
done

