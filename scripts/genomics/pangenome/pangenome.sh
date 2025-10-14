#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script implements pangenome analysis in two batches.
# The first one uses all strains, while the second one uses sinorhizobium only

# All strains
mkdir -p $folder_genomics/pangenome/all
for i in 2 3 4 5 6 8 9 10 11 13 15 16 17 19 20 {21..27} {29..37} {39..45}
do
    echo -e $folder_genomics/gff/g$i.gff
done >| $folder_genomics/pangenome/all/list_gffs.txt

mamba activate panaroo
panaroo \
    -i $folder_genomics/pangenome/all/list_gffs.txt \
    -o $folder_genomics/pangenome/all/ \
    -a core --aligner mafft \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes


# Sinorhizobium only
mkdir -p $folder_genomics/pangenome/sino
for i in 4 5 6 8 9 10 11 13 16 17 19 20 {21..27} {29..37} {39..45}
do
    echo -e $folder_genomics/gff/g$i.gff
done >| $folder_genomics/pangenome/sino/list_gffs.txt

mamba activate panaroo
panaroo \
    -i $folder_genomics/pangenome/sino/list_gffs.txt \
    -o $folder_genomics/pangenome/sino/ \
    -a core --aligner mafft \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes
