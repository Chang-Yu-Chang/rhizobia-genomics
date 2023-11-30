#!/usr/bin/env zsh
source ~/.zshrc

conda activate





"THIS IS NOT RIGHT. THIS FILE IS FOR ALIGNING raw reads to reference, not to find SNPs; create a new folder for calling that. Maybe 07-popgen_prep?"









# 6.2 Align Each Assembled Genome to the Reference
mamba activate minimap2
minimap2 -ax map-ont $1 "$folder_genomes/01-fasta/$2.fasta" > "$folder_genomes/$2/02-guided_assembly/.sam"

# 6.3 Cover sam to bam
mamba activate samtools
mkdir -p "$folder_anvio/06-alignment/bam"
# Convert SAM to BAM: For each alignment use samtools to convert the SAM file to BAM and sort it.
samtools view -bS "$folder_anvio/06-alignment/sam/$i.sam" | samtools sort -o "$folder_anvio/06-alignment/bam/$i.bam"
# Index the BAM file: Index the resulting BAM file
samtools index "$folder_anvio/06-alignment/bam/$i.bam"


# for i in Chang_Q5C_{1..19} usda1106 em1021 em1022 wsm419
# do
# done
