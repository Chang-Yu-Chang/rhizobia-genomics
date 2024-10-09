#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
mamba activate panaroo
mkdir -p $folder_genomics/pangenome/isolates

# 1. A total of 36 strains. Skip g38 and g40
# for i in {1..32} 34 36 37 38
# do
#     echo -e $folder_genomics/gff/$genome_ids[$i].gff
# done >| $folder_genomics/pangenome/isolates/list_gffs.txt
#
# panaroo \
#     -i $folder_genomics/pangenome/isolates/list_gffs.txt \
#     -o $folder_genomics/pangenome/isolates \
#     -a pan --aligner mafft \
#     --core_threshold 0.95 \
#     -t 10 \
#     --clean-mode strict \
#     --remove-invalid-genes

# 2. Elevation medicae. Total 10
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


# 3. Urbanization meliloti. Total 17
mkdir -p $folder_genomics/pangenome/urbn_mel
for i in {21..27} {31..37} 39 41 43
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



# cat $folder_genomics/pangenome/isolates/list_gffs.txt > $folder_genomics/pangenome/genomes/list_gffs.txt
#
# echo $folder_genomics/gff/genomes/em1021.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
# echo $folder_genomics/gff/genomes/em1022.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
# echo $folder_genomics/gff/genomes/usda1106.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
# echo $folder_genomics/gff/genomes/wsm419.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
# echo $folder_genomics/gff/genomes/casidaa.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
