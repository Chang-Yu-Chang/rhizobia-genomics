#!/usr/bin/env zsh
source ~/.zshrc

# This summarizes the pangenome result in a txt file
# $1: PAN.db
# $2: GENOMES.db
# $3: output directory
# $4: project name

conda activate
mamba activate anvio-8

# Add a default collection
anvi-script-add-default-collection \
    -p $1 \
    -c $2 \
    -C default \
    -b everything

anvi-summarize \
    --force-overwrite \
    -p $1 \
    -g $2 \
    -o $3 \
    -C default
# -C the collection of bins has to be specificed in the GUI using anvi-display-pan

# unzip the txt file
gunzip $3/$4'_gene_clusters_summary.txt.gz'
