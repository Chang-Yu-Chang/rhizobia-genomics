#!/usr/bin/env zsh
source ~/.zshrc

# This calls SVs via sniffle
# $1: target genome in bam
# $2: target genome in snf

mamba activate sniffles

sniffles --input $1 --snf $2 --allow-overwrite
