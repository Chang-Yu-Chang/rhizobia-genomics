#!/usr/bin/env zsh
source ~/.zshrc

# This annotates de novo assembled genomes via bakta
# $1: genome in fasta
# $2: bakta folder
# $3: bakta_database

conda activate
mamba activate bakta

bakta --db $3 --verbose --output $2 --thread 8 $1
# `--db "$bakta_database"` database directory
# `--verbose` Print verbose information
# `--output "$bakta_folder"` output directory
# `--thread 8`
# `"$medaka_consensus"`  Genome sequences in (zipped) fasta format


# Show the circos plot
bakta_plot --output $2 --prefix circos "$2/genome.json"
