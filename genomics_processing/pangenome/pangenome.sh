#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
mamba activate panaroo
mkdir -p $folder_genomics/pangenome/isolates
mkdir -p $folder_genomics/pangenome/genomes

# Create a list of isolate gff. Skip g28
for i in {1..22} {24..32}
do
    echo -e $folder_genomics/gff/genomes/$genome_ids[$i].gff
done >| $folder_genomics/pangenome/isolates/list_gffs.txt

panaroo \
    -i $folder_genomics/pangenome/isolates/list_gffs.txt \
    -o $folder_genomics/pangenome/isolates \
    -a core --aligner mafft \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes
# -t N_CPU, --threads N_CPU number of threads to use (default=1)


cat $folder_genomics/pangenome/isolates/list_gffs.txt > $folder_genomics/pangenome/genomes/list_gffs.txt

echo $folder_genomics/gff/genomes/em1021.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/genomes/em1022.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/genomes/usda1106.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/genomes/wsm419.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
echo $folder_genomics/gff/genomes/casidaa.gff >> $folder_genomics/pangenome/genomes/list_gffs.txt
