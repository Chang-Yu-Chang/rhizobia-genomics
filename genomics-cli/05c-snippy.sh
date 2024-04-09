#!/usr/bin/env zsh
source ~/.zshrc

# This calls SNPs via snippy
# $1: reference genome in fasta
# $2: raw reads in bam
# $3: snippy output folder

mamba activate snippy

snippy --cpus 10 --ref $1 --outdir $3 --force --bam $2
# --cpus N         Maximum number of CPU cores to use (default '8')
# --reference F    Reference genome. Supports FASTA, GenBank, EMBL (not GFF) (default '')
# --force          Force overwrite of existing output folder (default OFF)
# --bam F          Use this BAM file instead of aligning reads (default '')


