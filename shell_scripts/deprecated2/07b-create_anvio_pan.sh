#!/usr/bin/env zsh
source ~/.zshrc

# This prepares the anvio object `db` for pangenome
# $1: list of external genomes in txt
# $2: pan genome output

conda activate
mamba activate anvio-8

anvi-gen-genomes-storage \
    -e $1 \
    -o $2 \
    --gene-caller Prodigal
# -e
# --gene-caller uses the external gene caller included in each genome's db
# -o it must ends with -GENOMES.db
