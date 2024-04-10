#!/usr/bin/env zsh
source ~/.zshrc

# This uses panaroo to analyze pangenome
# $1: output panaroo folder
# $2: prokka annotation gffs

mamba activate panaroo

panaroo \
    -i $2 \
    -o $1 \
    -a core --aligner mafft \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes
# -t N_CPU, --threads N_CPU number of threads to use (default=1)