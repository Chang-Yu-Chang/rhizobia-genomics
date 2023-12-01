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
    --ref "$folder_genomics/genomes/usda1106/02-denovo_assembly/ncbi/genome.gbff" \
    Chang_Q5C_2 Chang_Q5C_3 Chang_Q5C_4 Chang_Q5C_5 Chang_Q5C_6 \
    Chang_Q5C_8 Chang_Q5C_9 Chang_Q5C_10 Chang_Q5C_11 Chang_Q5C_13 \
    Chang_Q5C_15 Chang_Q5C_16 Chang_Q5C_17 Chang_Q5C_19 \
    --prefix "$folder_genomics/popgen/usda1106_snippy/usda1106_snippy"

# Call structure variants
conda activate
mamba activate sniffles

cd "$folder_genomics/popgen/snippy_usda1106"
for i in {1..19}; do; echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.snf" >> list_snf.tsv; done
sniffles --input list_snf.tsv --vcf sv.vcf


## reference: wsm419. Merge the set of 10 Ensifer medicae genomes
conda activate
mamba activate snippy

mkdir -p "$folder_genomics/popgen/snippy_wsm419"
snippy-core --ref "$folder_genomics/genomes/wsm419/02-denovo_assembly/ncbi/genome.gbff" \
    Chang_Q5C_2 Chang_Q5C_3 Chang_Q5C_4 Chang_Q5C_5 Chang_Q5C_6 \
    Chang_Q5C_8 Chang_Q5C_9 Chang_Q5C_10 Chang_Q5C_11 Chang_Q5C_13 \
    Chang_Q5C_15 Chang_Q5C_16 Chang_Q5C_17 Chang_Q5C_19 \
    --prefix "$folder_genomics/popgen/wsm419_snippy/wsm419_snippy"
    # Chang_Q5C_4 Chang_Q5C_5 Chang_Q5C_6 Chang_Q5C_8 Chang_Q5C_9 \
    # Chang_Q5C_11 Chang_Q5C_13 Chang_Q5C_16 Chang_Q5C_17 Chang_Q5C_19 \


# Call structure variants
conda activate
mamba activate sniffles

cd "$folder_genomics/popgen/snippy_wsm419"
for i in {1..19}; do; echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.snf" >> list_snf.tsv; done
sniffles --input list_snf.tsv --vcf sv.vcf








