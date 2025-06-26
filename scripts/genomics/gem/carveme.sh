#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script

cd $folder_shell
mkdir -p $folder_genomics/gem/carveme/

for i in {1..38}
do
    echo $genome_ids[$i]
    mamba activate carveme
    carve \
        $folder_genomics/faa/$genome_ids[$i].faa \
        --output $folder_genomics/gem/carveme/$genome_ids[$i].xml \
        --gapfill LB \
        --init LB \
        --solver scip
done

