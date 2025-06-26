#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script

mkdir -p $folder_genomics/gem/simulate_lb/

for i in {1..38}
do
    echo $genome_ids[$i]
    mamba activate cobra
    python run_dfba.py \
        $folder_genomics/gem/carveme/$genome_ids[$i].xml \
        $folder_genomics/gem/simulate_lb/$genome_ids[$i].txt

done

for i in 11 12
do
    echo $genome_ids[$i]
    mamba activate cobra
    python remove_lb.py  --remove \
        $folder_genomics/gem/carveme/$genome_ids[$i].xml \
        $folder_genomics/gem/xml/$genome_ids[$i].xml
done
