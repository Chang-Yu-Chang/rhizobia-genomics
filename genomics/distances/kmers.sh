#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

cd $folder_shell

# This script creates and compares genomic kmer signatures

mkdir -p $folder_genomics/kmer/genomes

## Create kmer signatures
for i in {1..32}
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/fasta/genomes/$genome_ids[$i].fasta \
        $folder_genomics/kmer/genomes/$genome_ids[$i].sig
done

for ref in em1021 em1022 usda1106 wsm419 casidaa
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/fasta/genomes/$ref.fasta \
        $folder_genomics/kmer/genomes/$ref.sig
done

## Create a list of kmer signatures
for i in {1..32}
do
    echo $folder_genomics/kmer/genomes/$genome_ids[$i].sig
done |> $folder_genomics/kmer/genomes/list_sigs.txt

echo $folder_genomics/kmer/genomes/em1021.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/em1022.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/usda1106.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/wsm419.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/casidaa.sig >> $folder_genomics/kmer/genomes/list_sigs.txt


## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/kmer/genomes/list_sigs.txt \
    $folder_genomics/kmer/genomes/kmer.txt
