#!/usr/bin/env zsh
source ~/.zshrc
source ../00-env_vars.sh

# This script calls variants using snippy

cd $folder_shell
mkdir -p $folder_genomics/variants

for ref in em1021 wsm419; do;
    mkdir -p $folder_genomics/variants/$ref
    cd $folder_genomics/variants/$ref

    # Make list of tab
    for i in {1..32}; do;
        echo "$genome_ids[$i]\t$folder_genomics/fasta/genomes/$genome_ids[$i].fasta"
    done |> input.tab

    # Generate snippy scripts
    mamba activate snippy
    snippy-multi input.tab --ref $folder_genomics/fasta/genomes/$ref.fasta --cpus 16 > runme.sh

    # Run
    zsh runme.sh
done
