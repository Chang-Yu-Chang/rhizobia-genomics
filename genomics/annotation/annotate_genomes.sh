#!/usr/bin/env zsh
source ~/.zshrc
source ../00-env_vars.sh

# This script annotates genomes

cd $folder_shell

mkdir -p $folder_genomics/annotation/genomes

## Assembled genomes
for i in {1..32}
do
    echo $genome_ids[$i]
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/genomes/$genome_ids[$i]/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/genomes/$genome_ids[$i].fasta \
        $folder_genomics/annotation/genomes/$genome_ids[$i]/
done

## Reference genomes
for ref in em1021 em1022 usda1106 wsm419 casidaa
do
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$ref/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/genomes/$ref.fasta \
        $folder_genomics/annotation/genomes/$ref
done
