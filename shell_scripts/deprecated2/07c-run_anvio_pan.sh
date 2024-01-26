#!/usr/bin/env zsh
source ~/.zshrc

# This runs the pangenome analysis by anvio
# $1: pangenome. It ends with -GENOMES.db
# $2: project name
# $3: output directory
# $4: list of genomes included

conda activate
mamba activate anvio-8

anvi-pan-genome \
    --force-overwrite \
    -g $1 \
    -n $2 \
    --output-dir $3 \
    --num-threads 20 \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --genome-names $4 \
    --min-occurrence 1
# -g anvio genomes storage file
# -n project name
# --minbit the minimum minibit value
# --mcl-inflation MCL inflation parameter
# --genome-names use a subset of genomes
# --min-occurrence 2 removes singletons
# If you publish results from this workflow, please do not forget to cite DIAMOND (doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL (http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)
