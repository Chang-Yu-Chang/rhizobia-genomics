#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script consolidates genome fasta files into one folderx

mamba activate seqkit

for i in {1..38}
do
    dir=$folder_genomics/assembly/$genome_ids[$i]

    # Remove contigs < 100k
    seqkit seq -m 100000 -g $dir/medaka/consensus.fasta > $folder_genomics/fasta/genomes/$genome_ids[$i].fasta
done

