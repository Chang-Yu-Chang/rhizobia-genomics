#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script annotates the de novo assembled genomes

cd $folder_shell
#bakta_db=/Users/cychang/bioinformatics/bakta/db  # This database is mandatory and must be downloaded before annotation

mkdir -p $folder_genomics/annotation

for i in {1..38}
do
    genome_fa=$folder_genomics/genomes/$genome_ids[$i].fasta

    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$genome_ids[$i]/prokka
    zsh 05a-prokka.sh \
        $genome_fa \
        $folder_genomics/annotation/$genome_ids[$i]/prokka

    # # Annotate genomes via bakta
    # mkdir -p $folder_genomics/annotation/$genome_ids[$i]/batka
    # zsh 05b-bakta.sh \
    #     $genome_fa \
    #     $folder_genomics/annotation/$genome_ids[$i]/batka \
    #     $bakta_db

done
