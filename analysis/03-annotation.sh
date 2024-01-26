#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script annotates the de novo assembled genomes

cd $folder_shell

mkdir -p $folder_genomics/annotation

for i in {27..38}
do
    echo $genome_ids[$i]
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$genome_ids[$i]/prokka
    zsh 03a-prokka.sh \
        $folder_genomics/genomes/$genome_ids[$i].fasta \
        $folder_genomics/annotation/$genome_ids[$i]/prokka
done

for ref in em1021 em1022 usda1106 wsm419
do
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$ref/prokka
    zsh 03a-prokka.sh \
        $folder_genomics/genomes/$ref.fasta \
        $folder_genomics/annotation/$ref/prokka
done

# Consolidate annotation gff
mkdir -p $folder_genomics/gff
for i in {1..38}
do
    cp $folder_genomics/annotation/$genome_ids[$i]/prokka/annotated.gff $folder_genomics/gff/$genome_ids[$i].gff
done