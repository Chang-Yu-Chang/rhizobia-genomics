#!/usr/bin/env zsh
source ~/.zshrc

# This uses information acquired from the miniasm draft genome and re-downsamples the reads to ~100x coverage
# (do nothing if there isn't at least 100x coverage) with heavy weight applied to removing low quality reads
# (helps small plasmids stick around)
# $1: genome_size
# $2: draft_genome in fasta
# $3: filtered_reads in fastq.gz
# $4: downsampled_reads2 in fastq

conda activate
mamba activate filtlong

filtlong --target_bases $(($1 * 100)) --mean_q_weight 10 --assembly $2 $3 | > $4
# `--target_bases $((genome_size * 100))` allow the coverage X100 of genome size
# `--mean_q_weight 10` specify a mean quality weight of 10 (instead of the default of 1) makes mean read quality more important when choosing the best reads
# `--assembly` use reference assembly in FASTA format
