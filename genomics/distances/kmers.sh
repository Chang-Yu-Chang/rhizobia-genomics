#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

cd $folder_shell

# This script creates and compares genomic kmer signatures

mkdir -p $folder_genomics/distances/kmer/signatures

## Create kmer signatures
for i in {1..38}
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/fasta/genomes/$genome_ids[$i].fasta \
        $folder_genomics/distances/kmer/signatures/$genome_ids[$i].sig
done

## Create a list of kmer signatures
for i in {1..38}
do
    echo $folder_genomics/distances/kmer/signatures/$genome_ids[$i].sig
done |> $folder_genomics/distances/kmer/list_sigs.txt


## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/distances/kmer/list_sigs.txt \
    $folder_genomics/distances/kmer/kmer.txt
