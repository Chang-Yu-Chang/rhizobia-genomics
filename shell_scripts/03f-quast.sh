#!/usr/bin/env zshs
source ~/.zshrc

# This checks genome quality
# $1: medaka consensus
# $2: quast folder

conda activate
mamba activate quast

quast $1 -o $2
# `-o` output directory
