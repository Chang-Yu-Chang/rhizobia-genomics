#!/usr/bin/env zsh
source ~/.zshrc

# This imports the ANI data to the pangenome file
# $1: ani in txt
# $2: PAN.db

conda activate
mamba activate anvio-8


anvi-import-misc-data \
    $1 \
    -p $2 \
    --target-data-table layers \
    --target-data-group ANI_percent_identity
# -e A two-column TAB-delimited flat text file that lists anvi'o contigs databases. Same txt used to make the anvio database
# -p useful when a anivo pangenome is available. The ANI will be written into the pangenome database
# --program program for computing ANI
# Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, please do not forget to properly credit their work.
