#!/usr/bin/env zsh
source ~/.zshrc

# This checks genome quality
# $1: genome in fasta
# $2: prokka folder

conda activate
mamba activate prokka

prokka --force --outdir $2 --kingdom Bacteria --locustag $1 --prefix annotated --gcode 11 $1
# `--force` force overwriting existing output folder
# `--outdir` o utput folder
# `--kingdom`
# `--locustag` locus tax prefix
# `--prefix` filname output prefix
# `--gcode` genetic code / translation table (set if --kingdom is set)
