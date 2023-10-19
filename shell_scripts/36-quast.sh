cd
source ~/.zshrc

# check genome quality
mamba activate quast
mamba env list

bakta_consensus=$1
quast_folder=$2

quast $bakta_consensus -o $quast_folder
# `-o` output directory
