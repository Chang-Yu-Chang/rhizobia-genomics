#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script searches and downloads the NCBI genomes
mamba activate ncbi-datasets

# Download the 18 Ensifer genomes
table_file=$folder_data/raw/ensifer_ncbi.csv
list_accessions=("${(@f)$(cat $table_file | cut -d ',' -f 1)}")
list_species=("${(@f)$(cat $table_file | cut -d ',' -f 2)}")
list_strains=("${(@f)$(cat $table_file | cut -d ',' -f 3)}")

datasets download genome accession \
    $list_accessions \
    --include genome \
    --filename $folder_genomics/ensifer_ncbi.zip

# GCF_002197065.1 E. meliloti usda1106
# GCF_000006965.1 E. meliloti em1021
# GCF_013315775.1 E. meliloti em1022
# GCF_002197445.1 E. meliloti usda1021
# GCF_004004435.1 E. meliloti usda1022
# GCF_000017145.1 E. medicae wsm419
# GCF_021052565.1 E. medicae wsm1115
# GCF_025200695.1 E. medicae su277
# GCF_000697965.2 E. adhaerens casidaa
# GCF_009883655.1 E. adhaerens corn53
# GCF_020035495.1 E. adhaerens w2a
# GCF_020035515.1 E. adhaerens w2b
# GCF_013488225.1 E. mexicanum ittgr7
# GCF_017488845.2 E. canadensis t173
# GCF_000018545.1 E. fredii ngr234
# GCF_000219415.3 E. fredii gr64
# GCF_000705595.2 E. americanum ccgm7
# GCF_002288525.1 E. sojae ccbau05684
# GCF_008932245.1 E. alkalisoli yic4027

# Clean up the files and names
cd $folder_genomics
unzip ensifer_ncbi.zip
rm -rf ncbi_genomes
mv ncbi_dataset/data ncbi_genomes
rm -rf README.md ensifer_ncbi.zip ncbi_dataset

# Move the genome files
mkdir -p $folder_genomics/fasta/genomes
n_strains="${#list_strains}"
for i in {1..$n_strains}; do
    cp $folder_genomics/ncbi_genomes/$list_accessions[$i]/*_genomic.fna $folder_genomics/fasta/genomes/$list_strains[$i].fasta
done

# Remove the zip folder
rm -rf $folder_genomics/ncbi_genomes
