#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script aggregates the SNPs and SVs into one vcf

cd $folder_shell
echo "08-popgen"

# Snippy
## reference: usda1106. Merge the vcf of the 14 Ensifer genomes
#mkdir -p "$folder_anvio/06-alignment/snippy/core"
conda activate
mamba activate snippy

mkdir -p "$folder_genomics/popgen/snippy_usda1106"
snippy-core \
    --ref "$folder_genomes/usda1106/02-denovo_assembly/ncbi/genome.gbff" \
    "$folder_genomes/Chang_Q5C_2/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_3/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_4/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_5/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_6/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_8/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_9/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_10/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_11/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_13/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_15/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_16/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_17/03-variant_calling/snippy_usda1106/snippy" \
    "$folder_genomes/Chang_Q5C_19/03-variant_calling/snippy_usda1106/snippy" \
    --prefix "$folder_genomics/popgen/snippy_usda1106/snippy_usda1106"

# Call structure variants
conda activate
mamba activate sniffles

cd "$folder_genomics/popgen/snippy_usda1106"
for i in {2..6} {8..11} 13 {15..17} 19
do
    echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.snf"
done |> list_snf.tsv
sniffles --allow-overwrite --input list_snf.tsv --vcf sv.vcf


## reference: wsm419. Merge the set of 10 Ensifer medicae genomes
conda activate
mamba activate snippy

mkdir -p "$folder_genomics/popgen/snippy_wsm419"
snippy-core \
    --ref "$folder_genomes/usda1106/02-denovo_assembly/ncbi/genome.gbff" \
    "$folder_genomes/Chang_Q5C_2/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_3/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_4/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_5/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_6/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_8/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_9/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_10/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_11/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_13/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_15/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_16/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_17/03-variant_calling/snippy_wsm419/snippy" \
    "$folder_genomes/Chang_Q5C_19/03-variant_calling/snippy_wsm419/snippy" \
    --prefix "$folder_genomics/popgen/snippy_wsm419/snippy_wsm419"


# Call structure variants
conda activate
mamba activate sniffles

cd "$folder_genomics/popgen/snippy_wsm419"
for i in {2..6} {8..11} 13 {15..17} 19
do
    echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.snf"
done |> list_snf.tsv
sniffles --allow-overwrite --input list_snf.tsv --vcf sv.vcf







