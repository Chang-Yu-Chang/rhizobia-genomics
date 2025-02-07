#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script runs defense finder
mkdir -p $folder_genomics/defense/

for i in {1..38}
do
    echo $genome_ids[$i]
    mkdir -p $folder_genomics/defense/$genome_ids[$i]

    # Defense finder
    # mamba activate defense-finder
    # defense-finder run \
    #     $folder_genomics/faa/$genome_ids[$i].faa \
    #     -o $folder_genomics/defense/$genome_ids[$i]/defense_finder \
    #     -a

    # padloc
    mamba activate padloc
    mkdir -p $folder_genomics/defense/$genome_ids[$i]/padloc
    padloc \
        --faa $folder_genomics/faa/$genome_ids[$i].faa \
        --gff $folder_genomics/gff/$genome_ids[$i].gff \
        --outdir $folder_genomics/defense/$genome_ids[$i]/padloc \
        --cpu 4
done
