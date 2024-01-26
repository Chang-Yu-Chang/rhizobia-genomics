#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script annotates the de novo assembled genomes

cd $folder_shell

mkdir -p $folder_genomics/annotation

for i in 1 8 16 18 {22..38}
do
    genome_fa=$folder_genomics/genomes/$genome_ids[$i].fasta

    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$genome_ids[$i]/prokka
    zsh 03a-prokka.sh \
        $genome_fa \
        $folder_genomics/annotation/$genome_ids[$i]/prokka

done

# Consolidate annotation gff
mkdir -p $folder_genomics/gff
for i in {1..38}
do
    cp $folder_genomics/annotation/prokka/$genome_ids[$i]/annotation.gff $folder_genomics/gff/$genome_ids[$i].gff
done
