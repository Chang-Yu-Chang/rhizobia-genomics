#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script performs reference-guided genome assembly

cd $folder_shell

# # Download and prepared NCBI genomes
# zsh 03a-download_ncbi_genomes.sh \
#     "$folder_data/temp/00-genomes_mapping.csv" \
#     $folder_genomics
#
# # Index the reference genomes using the four NCBI genomes
# for i in {20..23}
# do
#     zsh 03b-index_references.sh \
#         "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
#         "$folder_genomics/references/$sample_ids[$i].mmi"
# done


# Process each genome
for i in {28..41}
do
    mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling"
    # Use reference E. meliloti usda1106
    ## Align the genome to the reference usda1106
    mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/snippy"
    zsh 03c-align_genome.sh \
        "$folder_genomics/references/usda1106.mmi" \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.sam"

    ## Convert SAM to BAM
    zsh 03d-convert_sam_to_bam.sh \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.sam" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.bam"

    # Call SNPs via snippy using reference usda1106
    zsh 03e-snippy.sh \
        "$folder_genomics/genomes/usda1106/02-denovo_assembly/ncbi/genome.gbff" \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/snippy"

    ## Call SVs via sniffle
    zsh 03f-sniffle.sh \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.bam" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_usda1106/genome.snf"



    # Use reference E. medicae wsm419
    ## Align the genome to the reference usda1106
    mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/snippy"
    zsh 03c-align_genome.sh \
        "$folder_genomics/references/usda1106.mmi" \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.sam"

    ## Convert SAM to BAM
    zsh 03d-convert_sam_to_bam.sh \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.sam" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.bam"

    # Call SNPs via snippy using reference usda1106
    zsh 03e-snippy.sh \
        "$folder_genomics/genomes/usda1106/02-denovo_assembly/ncbi/genome.gbff" \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/snippy"

    ## Call SVs via sniffle
    zsh 03f-sniffle.sh \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.bam" \
        "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_wsm419/genome.snf"



done

# Merge the called SNPs and SVs
# Call structure variants
#cd "$folder_anvio/06-alignment/sniffles"
# echo $folder_anvio/06-alignment/snf/$i.snf >> list_snf.tsv; done
# sniffles --input list_snf.tsv --vcf sv.vcf






