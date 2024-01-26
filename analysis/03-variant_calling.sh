#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script calls SNPs

cd $folder_shell
mkdir -p $folder_genomics/variants



# Process each genome
for ref in em1021 wsm419
do
    # Index the reference genome
    conda activate
    mamba activate minimap2
    minimap2 -d $folder_genomics/genomes/$ref.mmi $folder_genomics/genomes/$ref.fasta


    mkdir -p $folder_genomics/variants/$ref
    
    for i in {1..38}
    do
        echo $genome_ids[$i]
        mkdir -p $folder_genomics/variants/$ref/$genome_ids[$i]
        dir=$folder_genomics/variants/$ref/$genome_ids[$i]

        # Align the filtered reads to reference genome
        zsh 03a-align_reads.sh \
            $folder_genomics/genomes/$ref.mmi \
            $folder_genomics/raw_reads/$genome_ids[$i].fastq.gz \
            $dir/genome.sam

        # Convert SAM to BAM
        zsh 03b-convert_sam_to_bam.sh \
            $dir/genome.sam \
            $dir/genome.bam

        # Call SNPs via snippy using reference
        mkdir -p $dir/snippy
        zsh 03c-snippy.sh \
            $folder_genomics/genomes/$ref.fasta \
            $dir/genome.bam \
            $dir/snippy

        # zsh 03a-align_reads.sh \
        #     $folder_genomics/genomes/$ref.fasta \
        #     $folder_genomics/assembly/$genome_ids[$i]/filtered_reads.fastq.gz \
        #     $dir
    done

    # Snippy all samples
    snippy-core \
        --ref $folder_genomes/genomes/$ref.fasta \
        $folder_genomics/variants/* \
        --prefix $folder_genomics/variants_core/$ref/core
done


# # Call structure variants
# conda activate
# mamba activate sniffles

# mkdir -p "$folder_genomics/popgen/$ref/sniffle"
# cd "$folder_genomics/popgen/$ref/sniffle"
# for i in {2..6} {8..11} 13 {15..17} 19 {20..37}
# do
#     echo "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.snf"
# done |> list_snf.tsv
# sniffles --allow-overwrite --input list_snf.tsv --vcf sv.vcf


#         # Call structural variants
#         ## sniffe
#     done
# done


#    for ref in usda1106 wsm419
#     do
#         # Use reference E. meliloti usda1106  and E. medicae wsm419
#         ## Align the genome to the reference
#         mkdir -p "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/snippy"
#         zsh 03c-align_genome.sh \
#             "$folder_genomics/references/$ref.mmi" \
#             "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
#             "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.sam"

#         ## Convert SAM to BAM
#         zsh 03d-convert_sam_to_bam.sh \
#             "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.sam" \
#             "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.bam"

#         # Call SNPs via snippy using reference
#         zsh 03e-snippy.sh \
#             "$folder_genomics/genomes/$ref/02-denovo_assembly/genome.fasta" \
#             "$folder_genomics/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
#             "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/snippy"

#         ## Call SVs via sniffle
#         zsh 03f-sniffle.sh \
#             "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.bam" \
#             "$folder_genomics/genomes/$sample_ids[$i]/03-variant_calling/snippy_$ref/genome.snf"
#     done