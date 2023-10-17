# This script executes the bacterial genome assembly pipeline inspired by plasmisaurus
# https://www.plasmidsaurus.com/faq/#bactseq

# Run the scripts
cd ~/Desktop/lab/local-adaptation/analysis

zsh 31-filter_reads.sh
zsh 32-miniasm.sh
zsh 33-flye.sh
zsh 34-medaka.sh
zsh 35-bakta.sh
zsh 36-checkm.sh
zsh 37-mash.sh
zsh 38-sourmash.sh
