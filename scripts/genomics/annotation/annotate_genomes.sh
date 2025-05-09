#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script annotates genomes

cd $folder_shell

mkdir -p $folder_genomics/annotation

## Assembled genomes
for i in {1..38}
do
    echo $genome_ids[$i]
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$genome_ids[$i]/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/genomes/$genome_ids[$i].fasta \
        $folder_genomics/annotation/$genome_ids[$i]/
done

## Reference genomes
for ref in em1021 em1022 usda1106 wsm419 casidaa
#for ref in em1021 wsm419 ngr234
do
    # Annotate genomes via prokka
    mkdir -p $folder_genomics/annotation/$ref/
    zsh 03a-prokka.sh \
        $folder_genomics/fasta/ncbi/$ref.fasta \
        $folder_genomics/annotation/$ref
done
