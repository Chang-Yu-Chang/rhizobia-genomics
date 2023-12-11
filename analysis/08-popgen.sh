#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script aggregates the SNPs and SVs into one vcf

cd $folder_shell
echo "08-popgen"


## reference: usda1106. Merge the vcf of the 14 Ensifer genomes
for ref in usda1106 wsm419
do
    # Snippy
    conda activate
    mamba activate snippy

    mkdir -p "$folder_genomics/popgen/snippy_$ref"
    snippy-core \
        --ref "$folder_genomes/$ref/02-denovo_assembly/genome.fasta" \
        "$folder_genomes/Chang_Q5C_"{2..6}"/03-variant_calling/snippy_$ref/snippy" \
        "$folder_genomes/Chang_Q5C_"{8..11}"/03-variant_calling/snippy_$ref/snippy" \
        "$folder_genomes/Chang_Q5C_13/03-variant_calling/snippy_$ref/snippy" \
        "$folder_genomes/Chang_Q5C_"{15..17}"/03-variant_calling/snippy_$ref/snippy" \
        "$folder_genomes/Chang_Q5C_19/03-variant_calling/snippy_$ref/snippy" \
        "$folder_genomes/Chang_W8S_"{1..18}"/03-variant_calling/snippy_$ref/snippy" \
        --prefix "$folder_genomics/popgen/snippy_$ref/snippy_$ref"

    # Call structure variants
    conda activate
    mamba activate sniffles

    cd "$folder_genomics/popgen/snippy_$ref"
    for i in {2..6} {8..11} 13 {15..17} 19 {24..41}
    do
        echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.snf"
    done |> list_snf.tsv
    sniffles --allow-overwrite --input list_snf.tsv --vcf sv.vcf

done
