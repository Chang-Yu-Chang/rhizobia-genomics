#!/usr/bin/env zsh

folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
cd "$folder_data/temp/anvio"
ref_genome="$folder_data/temp/anvio/reference_genome.mmi"

# 1. Index the Reference Genome
# For long reads, tools like Minimap2 or NGMLR are often preferred
mamba activate minimap2
minimap2 -d $ref_genome "$folder_data/temp/anvio/contigs/em1021.fasta"

# 2, Align Each Assembled Genome to the Reference
mkdir -p "$folder_data/temp/anvio/sam"
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    minimap2 -ax map-ont $ref_genome "$folder_data/temp/anvio/contigs/$i.fasta" > "$folder_data/temp/anvio/sam/$i.sam"
done

# Cover sam to bam
mamba activate samtools
mkdir -p "$folder_data/temp/anvio/bam"
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    # Convert SAM to BAM: For each alignment use samtools to convert the SAM file to BAM and sort it.
    samtools view -bS "$folder_data/temp/anvio/sam/$i.sam" | samtools sort -o "$folder_data/temp/anvio/bam/$i.bam"
    # Index the BAM file: Index the resulting BAM file
    samtools index "$folder_data/temp/anvio/bam/$i.bam"
done

# i=Chang_Q5C_1
# samtools mpileup "$folder_data/temp/anvio/bam/$i.bam"

# 3. call SNPs variation
mamba activate anvio-8
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"

# Profile a BAM file. For EACH SAMPLE
mkdir -p "$folder_data/temp/anvio/single_profiles"
#i=Chang_Q5C_1
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    anvi-profile \
        -i "$folder_data/temp/anvio/bam/$i.bam" \
        -c "$folder_data/temp/anvio/genomes/$i.db" \
        -o "$folder_data/temp/anvio/profiles/$i" \
        --sample-name $i \
        --profile-SCVs --overwrite-output-destinations
done

# Merge all anvi’o profiles using the program anvi-merge
anvi-merge "$folder_data/temp/anvio/profiles/*/PROFILE.db" -o SAMPLES-MERGED -c contigs.db

anvi-gen-variability-profile \
    -p profile-db \
    -c contigs-db \
    -C DEFAULT \
    -b EVERYTHING







