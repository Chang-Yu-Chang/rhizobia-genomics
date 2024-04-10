#!/usr/bin/env zsh
source ~/.zshrc

# This calls variants from core genes
# $1: reference fasta
# $2: alignment
# $3: output call in bcf

# Make bam file
mamba activate samtools
samtools faidx core_gene_alignment.aln
samtools view -bt core_gene_alignment.aln.fai -o core_gene_alignment.bam core_gene_alignment.aln

# Sort and index the bam file
samtools sort -o core_gene_alignment_sorted.bam core_gene_alignment.bam
samtools index core_gene_alignment_sorted.bam

# Call variants using bcftools
mamba activate bcftools
bcftools mpileup -Ou -f core_gene.fasta core_gene_alignment_sorted.bam | bcftools call -mv -Ov -o core_gene_variants.vcf
