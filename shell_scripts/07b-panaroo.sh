#!/usr/bin/env zsh
source ~/.zshrc

# This uses panaroo to analyze pangenome
# $1: output panaroo folder
# $2: folder where the prokka annotation gffs are stored

conda activate
mamba activate panaroo

panaroo -i $2 \
    -o $1 \
    -t 10 \
    --clean-mode moderate --remove-invalid-genes
