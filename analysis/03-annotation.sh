#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script annotates the de novo assembled genomes

cd $folder_shell

mkdir -p $folder_genomics/annotation/genomes

# Annotate genomes
for i in {1..32}
do
    echo $genome_ids[$i]
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/genomes/$genome_ids[$i]/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/genomes/$genome_ids[$i].fasta \
        $folder_genomics/annotation/genomes/$genome_ids[$i]/
done

for ref in em1021 em1022 usda1106 wsm419 casidaa
do
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$ref/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/genomes/$ref.fasta \
        $folder_genomics/annotation/genomes/$ref
done

# Consolidate genome gff
mkdir -p $folder_genomics/gff/genomes
for i in {1..32}; do
    cp $folder_genomics/annotation/genomes/$genome_ids[$i]/annotated.gff $folder_genomics/gff/genomes/$genome_ids[$i].gff
done

for ref in em1021 em1022 usda1106 wsm419 casidaa; do
    cp $folder_genomics/annotation/genomes/$ref/annotated.gff $folder_genomics/gff/genomes/$ref.gff
done

