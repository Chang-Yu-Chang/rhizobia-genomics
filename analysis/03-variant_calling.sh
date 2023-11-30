#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script performs reference-guided genome assembly

cd $folder_shell
echo "03-variant calling"

# Download and prepared NCBI genomes
zsh 03a-download_ncbi_genomes.sh \
    "$folder_data/temp/00-genomes_mapping.csv" \
    $folder_genomics

# Index the reference genomes using the four NCBI genomes
for i in {20..23}
do
    zsh 02b-index_references.sh \
        "$folder_genomics/references/$sample_ids[$i].mmi" \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"
done



for i in {1..23}
do
    # Align the genome to the reference usda1106
    zsh 02c-align_genomes.sh \
        "$folder_genomics/references/usda1106.mmi" \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_anvio/06-alignment/sam/$i.sam"

    # Convert SAM to BAM
    zsh 02d-convert_sam_to_bam.sh

    # Call SNPs via snippy
    zsh 02e-snippy.sh

    # Call SVs via sniffle
    zsh 02f-sniffle.sh


    # zsh 02a-miniasm.sh  \
    #     "$folder_genomics/reference/usda1106.mmi" \
    #     $sample_id
    #
    #
    # zsh 02b-align_genomes.sh  \
    #     "$folder_genomics/reference/usda1106.mmi" \
    #     $sample_id

done


