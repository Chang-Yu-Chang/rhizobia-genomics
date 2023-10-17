# This script executes the bacterial genome assembly pipeline inspired by plasmisaurus
# https://www.plasmidsaurus.com/faq/#bactseq

# Run the scripts
zsh 31-filter_reads.sh
zsh 32-miniasm.sh
zsh 33-flye_medaka.sh
