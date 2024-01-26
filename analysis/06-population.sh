#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

cd $folder_shell

# This script performs analysis comparing distance between genomes
mkdir -p $folder_genomics/popgen


# 1. Calculate ANI
mkdir -p $folder_genomics/popgen/fastani
for i in {1..42}
do
    echo -e $folder_genomics/genomes/$genome_ids[$i].fasta
done >| $folder_genomics/popgen/fastani/list_genomes_for_ani.txt

zsh 06a-fastani.sh \
    $folder_genomics/popgen/fastani/list_genomes_for_ani.txt \
    $folder_genomics/popgen/fastani/ani.txt

# 2. Compare the k-mer signature among genomes
mkdir -p $folder_genomics/popgen/kmer

## Create kmer signatures
for i in {1..38}
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/genomes/$genome_ids[$i].fasta \
        $folder_genomics/popgen/kmer/$genome_ids[$i].sig
done

for ref in em1021 em1022 usda1106 wsm419
do
    zsh 06b-genome_kmer.sh \
        $folder_genomics/genomes/$ref.fasta \
        $folder_genomics/popgen/kmer/$ref.sig
done

## Create a list of kmer signatures
for i in {1..38}; do; echo $folder_genomics/popgen/kmer/$genome_ids[$i].sig
done |> $folder_genomics/popgen/kmer/list_sig.txt
echo $folder_genomics/popgen/kmer/em1021.sig >> $folder_genomics/popgen/kmer/list_sig.txt
echo $folder_genomics/popgen/kmer/em1022.sig >> $folder_genomics/popgen/kmer/list_sig.txt
echo $folder_genomics/popgen/kmer/usda1106.sig >> $folder_genomics/popgen/kmer/list_sig.txt
echo $folder_genomics/popgen/kmer/wsm419.sig >> $folder_genomics/popgen/kmer/list_sig.txt


## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/popgen/kmer/list_sig.txt \
    $folder_genomics/popgen/kmer/kmer.txt
    

# 3. aggregates the SNPs and SVs into one vcf. This is using the raw read alignment
for ref in em1021 wsm419
do
    mkdir -p $folder_genomics/popgen/snp/$ref
    # Snippy to call SNPs
    mamba activate snippy
    mkdir -p $folder_genomics/popgen/snp/$ref/
    # for i in {4..6} {8..11} 13 {16..17} 19 {20..37}
    # do
    #     mkdir -p $folder_genomics/popgen/read_$ref/snippy/genomes/g$i
    #     cp $folder_genomes/$sample_ids[$i]/03-variant_calling/read_$ref/snippy/* $folder_genomics/popgen/read_$ref/snippy/genomes/g$i/
    # done
    snippy-core \
        --ref $folder_genomics/genomes/$ref.fasta \
        $folder_genomics/variants/$ref/em1021/* \
        --prefix  $folder_genomics/popgen/snp/$ref/core
done

# ref=em1021
# cd $folder_genomics/popgen/read_$ref/snippy/genomes/
# snippy-core \
#     --ref "$folder_genomes/$ref/02-denovo_assembly/genome.fasta" \
#     g20 g21 g22 g23 g24 g25 g26 g27 g31 g32 g33 g34 g35 g36 g37 \
#     --prefix "$folder_genomics/popgen/read_$ref/snippy/core"

# ref=wsm419
# cd $folder_genomics/popgen/read_$ref/snippy/genomes/
# snippy-core \
#     --ref "$folder_genomes/$ref/02-denovo_assembly/genome.fasta" \
#     g4 g5 g6 g8 g9 g11 g13 g16 g17 g19 g28 g29 g30 \
#     --prefix "$folder_genomics/popgen/read_$ref/snippy/core"






# # 3. aggregates the SNPs and SVs into one vcf. THis is using the genome alignment
# for ref in usda1106 wsm419
# do
#     mkdir -p "$folder_genomics/popgen/$ref"
#     # Snippy to call SNPs
#     conda activate
#     mamba activate snippy

#     mkdir "$folder_genomics/popgen/$ref/snippy/genomes"
#     for i in {2..6} {8..11} 13 {15..17} 19 {20..37}
#     do
#         mkdir -p $folder_genomics/popgen/$ref/snippy/genomes/g$i
#         cp $folder_genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/snippy/* $folder_genomics/popgen/$ref/snippy/genomes/g$i/
#     done


#     snippy-core \
#         --ref "$folder_genomes/$ref/02-denovo_assembly/genome.fasta" \
#         $folder_genomics/popgen/$ref/snippy/genomes/* \
#         --prefix "$folder_genomics/popgen/$ref/snippy/core"

#         # "$folder_genomes/Chang_Q5C_"{2..6z}"/03-variant_calling/snippy_$ref/snippy" \
#         # "$folder_genomes/Chang_Q5C_"{8..11}"/03-variant_calling/snippy_$ref/snippy" \
#         # "$folder_genomes/Chang_Q5C_13/03-variant_calling/snippy_$ref/snippy" \
#         # "$folder_genomes/Chang_Q5C_"{15..17}"/03-variant_calling/snippy_$ref/snippy" \
#         # "$folder_genomes/Chang_Q5C_19/03-variant_calling/snippy_$ref/snippy" \
#         # "$folder_genomes/Chang_W8S_"{1..18}"/03-variant_calling/snippy_$ref/snippy" \

#     # Call structure variants
#     conda activate
#     mamba activate sniffles

#     mkdir -p "$folder_genomics/popgen/$ref/sniffle"
#     cd "$folder_genomics/popgen/$ref/sniffle"
#     for i in {2..6} {8..11} 13 {15..17} 19 {20..37}
#     do
#         echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.snf"
#     done |> list_snf.tsv
#     sniffles --allow-overwrite --input list_snf.tsv --vcf sv.vcf


#     # # Filtering with bcftools
#     # conda activate
#     # mamba activate bcftools
#     #
#     # bcftools stats "$folder_genomics/popgen/$ref/snippy/core.vcf"
#     #
#     # bcftools filter -i'%QUAL>=30' input.vcf -o output_filtered.vcf

#     # LD pruning
#     # conda activate
#     # mamba activate vcftools
#     #
#     # mkdir -p "$folder_genomics/popgen/$ref/pruning"
#     # cd "$folder_genomics/popgen/$ref/pruning"
#     # vcftools --plink \
#     #     --vcf "$folder_genomics/popgen/$ref/snippy/core.vcf" \
#     #     --out "$folder_genomics/popgen/$ref/pruning/prepruned"

#     # # Convert PLINK files to binary format
#     # conda activate
#     # mamba activate plink2
#     #
#     # plink2 \
#     #     --vcf "$folder_genomics/popgen/$ref/snippy/core.vcf" \
#     #     --make-bed --allow-extra-chr \
#     #     --out "$folder_genomics/popgen/$ref/pruning/prepruned"
#     #
#     # # Perform LD pruning
#     # # plink2 \
#     # #     --bfile "$folder_genomics/popgen/$ref/pruning/prepruned" \
#     # #     --allow-extra-chr --bad-ld \
#     # #     --indep-pairwise 1000 100 0.2 \
#     # #     --out "$folder_genomics/popgen/$ref/pruning/pruned.prune.in"
#     #
#     # # Extract LD-pruned variants from VCF
#     # conda activate
#     # mamba activate vcftools
#     #
#     # vcftools --recode \
#     #     --vcf "$folder_genomics/popgen/$ref/snippy/core.vcf" \
#     #     --positions "$folder_genomics/popgen/$ref/pruning/pruned.prune.in" \
#     #     --out "$folder_genomics/popgen/$ref/pruning/pruned.vcf"



#     # # Convert the vcf into a table
#     # conda activate
#     # mamba activate bcftools
#     # bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t]\n' your_file.vcf > output_table.tsv

# done


# # 4. Filter SNPs to retain those with delpth values 20-230
# mamba activate bcftools
# bcftools stats core.vcf

# mamba activate vcftools
# vcftools --vcf core.vcf --out output


















