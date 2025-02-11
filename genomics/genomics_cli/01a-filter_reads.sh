#!/usr/bin/env zsh
source ~/.zshrc

# This throws away the worst 5% reads
# $1: raw_reads
# $2: filtered_reads

mamba activate filtlong

filtlong --keep_percent 95 $1 | gzip > $2
#filtlong --keep_percent 95 $1 
# `--keep_percent 95` throw out the worst 5% of reads
