#!/usr/bin/env zshs
source ~/.zshrc

# This checks genome quality
# $1: medaka consensus
# $2: busco folder

conda activate
mamba activate busco

busco -i $1 --out_path $2 -o "results" -m genome --auto-lineage-prok
# `-i` input sequence file in FASTA format
# `-out_path` output folder name
# `-o` output project name
# `-m genome` BUSCO analysis mode. We use genome mode
# `--auto-lineage-prok` Automated lineage selection



