#!/usr/bin/env zsh
source ~/.zshrc
source ../00-env_vars.sh

cd $folder_shell

# This script performs analysis comparing distance between genomes/contigs

# 1. Calculate genomic ANI
mkdir -p $folder_genomics/fastani

## Create a list of genome fasta files
for i in {1..38}; do; echo -e $folder_genomics/genomes/$genome_ids[$i].fasta
done >| $folder_genomics/fastani/list_genomes.txt
echo $folder_genomics/genomes/em1021.fasta >> $folder_genomics/fastani/list_genomes.txt
echo $folder_genomics/genomes/em1022.fasta >> $folder_genomics/fastani/list_genomes.txt
echo $folder_genomics/genomes/usda1106.fasta >> $folder_genomics/fastani/list_genomes.txt
echo $folder_genomics/genomes/wsm419.fasta >> $folder_genomics/fastani/list_genomes.txt

## Compute ani
zsh 06a-fastani.sh \
    $folder_genomics/fastani/list_genomes.txt \
    $folder_genomics/fastani/ani.txt


# 2. Compare genomic kmer
mkdir -p $folder_genomics/kmer/genomes

## Create kmer signatures
for i in {1..38}
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/genomes/$genome_ids[$i].fasta \
        $folder_genomics/kmer/genomes/$genome_ids[$i].sig
done

for ref in em1021 em1022 usda1106 wsm419
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/genomes/$ref.fasta \
        $folder_genomics/kmer/genomes/$ref.sig
done

## Create a list of kmer signatures
for i in {1..38}; do; echo $folder_genomics/kmer/genomes/$genome_ids[$i].sig
done |> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/em1021.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/em1022.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/usda1106.sig >> $folder_genomics/kmer/genomes/list_sigs.txt
echo $folder_genomics/kmer/genomes/wsm419.sig >> $folder_genomics/kmer/genomes/list_sigs.txt


## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/kmer/genomes/list_sigs.txt \
    $folder_genomics/kmer/genomes/kmer.txt

# 3. Compute contig ANI (>10kb)


# 4. Compare contig kmer (>10kb)
mkdir -p $folder_genomics/kmer/contigs
list_contigs=($(ls $folder_genomics/contigs |  sed 's/\.fasta$//'))

## Create kmer signatures
for con in "${list_contigs[@]}";
do
    echo $con
    zsh 06b-genome_kmer.sh \
        $folder_genomics/contigs/$con.fasta \
        $folder_genomics/kmer/contigs/$con.sig
done

## Create a list of kmer signatures
list_sigs=($(ls $folder_genomics/kmer/contigs |  sed 's/\.sig$//'))
for sig in ${list_sigs[@]}; do;
    echo $folder_genomics/kmer/contigs/$sig.sig
done |> $folder_genomics/kmer/contigs/list_sigs.txt

echo $folder_genomics/popgen/kmer/em1021.sig >> $folder_genomics/kmer/contigs/list_sigs.txt
echo $folder_genomics/popgen/kmer/em1022.sig >> $folder_genomics/kmer/contigs/list_sigs.txt
echo $folder_genomics/popgen/kmer/usda1106.sig >> $folder_genomics/kmer/contigs/list_sigs.txt
echo $folder_genomics/popgen/kmer/wsm419.sig >> $folder_genomics/kmer/contigs/list_sigs.txt

## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/kmer/contigs//list_sigs.txt \
    $folder_genomics/kmer/contigs//kmer.txt

