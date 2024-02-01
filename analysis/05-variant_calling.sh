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
    mamba activate minimap2
    minimap2 -d $folder_genomics/genomes/$ref.mmi $folder_genomics/genomes/$ref.fasta
    
    mkdir -p $folder_genomics/variants/$ref
    
    for i in {1..38}
    do
        echo $genome_ids[$i]
        mkdir -p $folder_genomics/variants/$ref/$genome_ids[$i]
        dir=$folder_genomics/variants/$ref/$genome_ids[$i]

        # Align the filtered reads to reference genome
        zsh 05a-align_reads.sh \
            $folder_genomics/genomes/$ref.mmi \
            $folder_genomics/raw_reads/$genome_ids[$i].fastq.gz \
            $dir/genome.sam

        # Convert SAM to BAM
        zsh 05b-convert_sam_to_bam.sh \
            $dir/genome.sam \
            $dir/genome.bam

        # Call SNPs via snippy using reference
        mkdir -p $dir/snippy
        zsh 05c-snippy.sh \
            $folder_genomics/genomes/$ref.fasta \
            $dir/genome.bam \
            $dir/snippy
    done
done
