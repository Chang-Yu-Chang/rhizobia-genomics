#!/usr/bin/env zsh
source ~/.zshrc

# This converts sam to bam file
# $1: target genome sam file
# $2: target genome bam file

mamba activate samtools

# Convert SAM to BAM: For each alignment use samtools to convert the SAM file to BAM and sort it.
samtools view -bS $1| samtools sort -o $2
# Index the BAM file: Index the resulting BAM file
samtools index $2

