#!/usr/bin/env zsh
source ~/.zshrc

# This extract the filtered reads into a txt file
# $1: filtered_reads in fastq.gz
# $2: filtered_reads in txt

mamba activate bioawk

bioawk -c fastx '{print $name, $qual, length($seq)}' $1 > $2
