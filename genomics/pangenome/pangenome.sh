#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
mamba activate panaroo

# Elevation medicae. Total 10
mkdir -p $folder_genomics/pangenome/elev_med
for i in 4 5 6 8 9 11 13 16 17 19
do
    echo -e $folder_genomics/gff/g$i.gff
done >| $folder_genomics/pangenome/elev_med/list_gffs.txt

panaroo \
    -i $folder_genomics/pangenome/elev_med/list_gffs.txt \
    -o $folder_genomics/pangenome/elev_med \
    -a core --aligner mafft \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes

# Output updated GFF files
panaroo-generate-gffs \
    -i $folder_genomics/pangenome/elev_med/list_gffs.txt \
    -o $folder_genomics/pangenome/elev_med


# Urbanization meliloti Total 20
mkdir -p $folder_genomics/pangenome/urbn_mel
for i in {21..27} {31..37} 39 41 42 43 44 45
do
    echo -e $folder_genomics/gff/g$i.gff
done >| $folder_genomics/pangenome/urbn_mel/list_gffs.txt

panaroo \
    -i $folder_genomics/pangenome/urbn_mel/list_gffs.txt \
    -o $folder_genomics/pangenome/urbn_mel \
    -a core --aligner mafft \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes

panaroo-generate-gffs \
    -i $folder_genomics/pangenome/urbn_mel/list_gffs.txt \
    -o $folder_genomics/pangenome/urbn_mel

