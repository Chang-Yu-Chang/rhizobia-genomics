#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script performs reference-guided genome assembly

cd $folder_shell

# Download and prepared NCBI genomes
# zsh 03a-download_ncbi_genomes.sh \
#     "$folder_data/temp/00-genomes_mapping.csv" \
#     $folder_genomics

#Index the reference genomes using the four NCBI genomes
for i in {38..41}
do
    cp $folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta $folder_genomics/references/$sample_ids[$i].fasta
    zsh 03b-index_references.sh \
        "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        "$folder_genomics/references/$sample_ids[$i].mmi"
done


# Process each genome
#for i in {1..41}
for i in {2..6} {8..11} 13 {15..17} {19..41}
do
    mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling"

    # for ref in usda1106 wsm419
    # do
    #     # Use reference E. meliloti usda1106  and E. medicae wsm419
    #     ## Align the genome to the reference
    #     mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/snippy"
    #     zsh 03c-align_genome.sh \
    #         "$folder_genomics/references/$ref.mmi" \
    #         "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
    #         "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.sam"

    #     ## Convert SAM to BAM
    #     zsh 03d-convert_sam_to_bam.sh \
    #         "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.sam" \
    #         "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.bam"

    #     # Call SNPs via snippy using reference
    #     zsh 03e-snippy.sh \
    #         "$folder_genomics/genomes/$ref/02-denovo_assembly/genome.fasta" \
    #         "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
    #         "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/snippy"

    #     ## Call SVs via sniffle
    #     zsh 03f-sniffle.sh \
    #         "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.bam" \
    #         "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.snf"
    # done
    
    for ref in em1021 usda1106 wsm419
    do
        mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/read_$ref/snippy"
        ## Align the filtered reads to reference genomes
        zsh 03g-align_reads.sh \
            "$folder_genomics/references/$ref.fasta" \
            "$folder_genomics/genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.fastq.gz" \
            "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/read_$ref/genome.sam"
        
        ## Convert SAM to BAM
        zsh 03d-convert_sam_to_bam.sh \
            "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/read_$ref/genome.sam" \
            "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/read_$ref/genome.bam"

        # Call SNPs via snippy using reference
        zsh 03e-snippy.sh \
            "$folder_genomics/genomes/$ref/02-denovo_assembly/genome.fasta" \
            "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
            "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/read_$ref/snippy"
    done

done
