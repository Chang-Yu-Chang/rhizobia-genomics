#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script aggregates the SNPs and SVs into one vcf

cd $folder_shell
echo "08-popgen"


for ref in usda1106 wsm419
do
    mkdir -p "$folder_genomics/popgen/$ref"
    # Snippy to call SNPs
    conda activate
    mamba activate snippy

    mkdir "$folder_genomics/popgen/$ref/snippy/genomes"
    for i in {2..6} {8..11} 13 {15..17} 19 {24..41}
    do
        mkdir -p $folder_genomics/popgen/$ref/snippy/genomes/g$i
        cp $folder_genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/snippy/* $folder_genomics/popgen/$ref/snippy/genomes/g$i/
    done


    snippy-core \
        --ref "$folder_genomes/$ref/02-denovo_assembly/genome.fasta" \
        $folder_genomes/Chang_Q5C_"{2..6z}"/03-variant_calling/snippy_$ref/ \
        --prefix "$folder_genomics/popgen/$ref/snippy/core"

        # "$folder_genomes/Chang_Q5C_"{2..6z}"/03-variant_calling/snippy_$ref/snippy" \
        # "$folder_genomes/Chang_Q5C_"{8..11}"/03-variant_calling/snippy_$ref/snippy" \
        # "$folder_genomes/Chang_Q5C_13/03-variant_calling/snippy_$ref/snippy" \
        # "$folder_genomes/Chang_Q5C_"{15..17}"/03-variant_calling/snippy_$ref/snippy" \
        # "$folder_genomes/Chang_Q5C_19/03-variant_calling/snippy_$ref/snippy" \
        # "$folder_genomes/Chang_W8S_"{1..18}"/03-variant_calling/snippy_$ref/snippy" \
    # Call structure variants
    conda activate
    mamba activate sniffles

    mkdir -p "$folder_genomics/popgen/$ref/sniffle"
    cd "$folder_genomics/popgen/$ref/sniffle"
    for i in {2..6} {8..11} 13 {15..17} 19 {24..41}
    do
        echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.snf"
    done |> list_snf.tsv
    sniffles --allow-overwrite --input list_snf.tsv --vcf sv.vcf


    # Filtering with bcftools
    conda activate
    mamba activate bcftools

    bcftools stats "$folder_genomics/popgen/$ref/snippy/core.vcf"

    bcftools filter -i'%QUAL>=30' input.vcf -o output_filtered.vcf

    # LD pruning
    # conda activate
    # mamba activate vcftools
    #
    # mkdir -p "$folder_genomics/popgen/$ref/pruning"
    # cd "$folder_genomics/popgen/$ref/pruning"
    # vcftools --plink \
    #     --vcf "$folder_genomics/popgen/$ref/snippy/core.vcf" \
    #     --out "$folder_genomics/popgen/$ref/pruning/prepruned"

    # # Convert PLINK files to binary format
    # conda activate
    # mamba activate plink2
    #
    # plink2 \
    #     --vcf "$folder_genomics/popgen/$ref/snippy/core.vcf" \
    #     --make-bed --allow-extra-chr \
    #     --out "$folder_genomics/popgen/$ref/pruning/prepruned"
    #
    # # Perform LD pruning
    # # plink2 \
    # #     --bfile "$folder_genomics/popgen/$ref/pruning/prepruned" \
    # #     --allow-extra-chr --bad-ld \
    # #     --indep-pairwise 1000 100 0.2 \
    # #     --out "$folder_genomics/popgen/$ref/pruning/pruned.prune.in"
    #
    # # Extract LD-pruned variants from VCF
    # conda activate
    # mamba activate vcftools
    #
    # vcftools --recode \
    #     --vcf "$folder_genomics/popgen/$ref/snippy/core.vcf" \
    #     --positions "$folder_genomics/popgen/$ref/pruning/pruned.prune.in" \
    #     --out "$folder_genomics/popgen/$ref/pruning/pruned.vcf"



    # # Convert the vcf into a table
    # conda activate
    # mamba activate bcftools
    # bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t]\n' your_file.vcf > output_table.tsv

done

