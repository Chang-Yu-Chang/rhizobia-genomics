#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

cd $folder_shell

# This script performs analysis comparing distance between genomes
mkdir -p $folder_genomics/popgen

# 1. Calculate ANI
mkdir -p $folder_genomics/popgen/fastani

## Create a list of genome fasta files
for i in {1..38}; do; echo -e $folder_genomics/genomes/$genome_ids[$i].fasta
done >| $folder_genomics/popgen/fastani/list_genomes.txt
echo $folder_genomics/genomes/em1021.fasta >> $folder_genomics/popgen/fastani/list_genomes.txt
echo $folder_genomics/genomes/em1022.fasta >> $folder_genomics/popgen/fastani/list_genomes.txt
echo $folder_genomics/genomes/usda1106.fasta >> $folder_genomics/popgen/fastani/list_genomes.txt
echo $folder_genomics/genomes/wsm419.fasta >> $folder_genomics/popgen/fastani/list_genomes.txt

## Compute ani
zsh 06a-fastani.sh \
    $folder_genomics/popgen/fastani/list_genomes.txt \
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
done |> $folder_genomics/popgen/kmer/list_sigs.txt
echo $folder_genomics/popgen/kmer/em1021.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
echo $folder_genomics/popgen/kmer/em1022.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
echo $folder_genomics/popgen/kmer/usda1106.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
echo $folder_genomics/popgen/kmer/wsm419.sig >> $folder_genomics/popgen/kmer/list_sigs.txt


## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/popgen/kmer/list_sigs.txt \
    $folder_genomics/popgen/kmer/kmer.txt
    
# 3. aggregae the VCFs into one 
mkdir -p $folder_genomics/variants/vcfgz_wsm419
for i in {1..38}; do;
    bgzip -c $folder_genomics/variants/wsm419/$genome_ids[$i]/snippy/snps.vcf > $folder_genomics/variants/vcfgz_wsm419/$genome_ids[$i].vcf.gz
    tabix -p vcf $folder_genomics/variants/vcfgz_wsm419/$genome_ids[$i].vcf.gz
done

for i in {1..38}; do; echo $folder_genomics/variants/vcfgz_wsm419/$genome_ids[$i].vcf.gz
done |> $folder_genomics/variants/list_vcfgz_wsm419.txt

mamba activate bcftools
bcftools merge -O v --file-list $folder_genomics/variants/list_vcfgz_wsm419.txt \
    --force-samples \
    -o $folder_genomics/variants/wsm419.vcf.gz
    #$folder_genomics/variants/list_vcfs_wsm419.txt

bcftools merge -O v -o $folder_genomics/variants/wsm419.vcf -l $folder_genomics/variants/list_vcfs_wsm419.txt


# 4. Compare ani across contigs (>10kb)
# 5. Compare kmer across contigs (>10kb)
mkdir -p $folder_genomics/popgen/kmer_contigs

list_contigs=($(ls $folder_genomics/contigs |  sed 's/\.fasta$//'))

## Create kmer signatures
for con in "${list_contigs[@]}";
do
    echo $con
    zsh 06b-genome_kmer.sh \
        $folder_genomics/contigs/$con.fasta \
        $folder_genomics/popgen/kmer_contigs/$con.sig
done

## Create a list of kmer signatures
list_sigs=($(ls $folder_genomics/popgen/kmer_contigs |  sed 's/\.sig$//'))
for sig in ${list_sigs[@]}; do; 
    echo $folder_genomics/popgen/kmer_contigs/$sig.sig
done |> $folder_genomics/popgen/kmer_contigs/list_sigs.txt
# echo $folder_genomics/popgen/kmer/em1021.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
# echo $folder_genomics/popgen/kmer/em1022.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
# echo $folder_genomics/popgen/kmer/usda1106.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
# echo $folder_genomics/popgen/kmer/wsm419.sig >> $folder_genomics/popgen/kmer/list_sigs.txt
## Compare signatures
zsh 06c-compare_kmer.sh \
    $folder_genomics/popgen/kmer_contigs/list_sigs.txt \
    $folder_genomics/popgen/kmer_contigs/kmer_contigs.txt
    



# # 3. aggregates the SNPs and SVs into one vcf. This is using the raw read alignment
# for ref in em1021 wsm419
# do
#     mkdir -p $folder_genomics/popgen/snp/$ref
#     # Snippy to call SNPs
#     mamba activate snippy
#     mkdir -p $folder_genomics/popgen/snp/$ref/
#     # for i in {4..6} {8..11} 13 {16..17} 19 {20..37}
#     # do
#     #     mkdir -p $folder_genomics/popgen/read_$ref/snippy/genomes/g$i
#     #     cp $folder_genomes/$sample_ids[$i]/03-variant_calling/read_$ref/snippy/* $folder_genomics/popgen/read_$ref/snippy/genomes/g$i/
#     # done
#     snippy-core \
#         --ref $folder_genomics/genomes/$ref.fasta \
#         $folder_genomics/variants/$ref/*/snippy \
#         --prefix  $folder_genomics/popgen/snp/$ref/core
# done

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











