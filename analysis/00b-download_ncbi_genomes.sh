#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script searches and downloads the NCBI genomes
mamba activate ncbi-datasets

# Download the Ensifer genomes
datasets download genome accession \
    GCF_002197065.1 GCF_000006965.1 GCF_013315775.1 GCF_000017145.1 \
    --include genome \
    --filename $folder_genomics/ensifer_ncbi.zip
# GCF_002197065.1 is E. meliloti usda1106
# GCF_000006965.1 is E. meliloti em1021
# GCF_013315775.1 is E. meliloti em1022
# GCF_000017145.1 is E. medicae wsm419

# Clean up the files and names
cd $folder_genomics
unzip ensifer_ncbi.zip
rm -rf ncbi_genomes
mv ncbi_dataset/data ncbi_genomes
rm -rf README.md ensifer_ncbi.zip ncbi_dataset

# Move the genome files
mkdir -p $folder_genomics/genomes
cp $folder_genomics/ncbi_genomes/GCF_002197065.1/*_genomic.fna $folder_genomics/genomes/usda1106.fasta
cp $folder_genomics/ncbi_genomes/GCF_000006965.1/*_genomic.fna $folder_genomics/genomes/em1021.fasta
cp $folder_genomics/ncbi_genomes/GCF_013315775.1/*_genomic.fna $folder_genomics/genomes/em1022.fasta
cp $folder_genomics/ncbi_genomes/GCF_000017145.1/*_genomic.fna $folder_genomics/genomes/wsm419.fasta

# Remove the zip folder
rm -rf $folder_genomics/ncbi_genomes
