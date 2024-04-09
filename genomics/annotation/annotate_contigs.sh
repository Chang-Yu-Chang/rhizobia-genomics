#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script annotates contigs

cd $folder_shell

mkdir -p $folder_genomics/annotation/contigs

# Annotate contigs
contigs=($(ls $folder_genomics/fasta/contigs/ | sed 's/\.fasta//'))
echo $contigs

for i in {1..${#contigs[@]}}
#for i in {1..9};
do
    echo $contigs[$i]
    # Annotate contigs via prokka
    mkdir -p $folder_genomics/annotation/contigs/$contigs[$i]/

    zsh 03a-prokka.sh \
        $folder_genomics/fasta/contigs/$contigs[$i].fasta \
        $folder_genomics/annotation/contigs/$contigs[$i]
done


# Consolidate contig gff
mkdir -p $folder_genomics/gff/contigs
for i in {1..${#contigs[@]}};
do
    cp $folder_genomics/annotation/contigs/$genome_ids[$i]/annotated.gff $folder_genomics/gff/contigs/$genome_ids[$i].gff
done
