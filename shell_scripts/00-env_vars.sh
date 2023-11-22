#!/usr/bin/env zshs

# This scripts save the environment varaibles shared by all shell scripts

source ~/.zshrc
conda activate

# Path names
folder_shell="/Users/cychang/Desktop/lab/local-adaptation/" # # Change this to your project folder
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data" # Change this to your data folder
folder_raw="$folder_data/raw/"
folder_temp="$folder_data/temp/"
folder_genomics="$folder_data/genomics/"
folder_genomes="$folder_data/genomics/genomes"

mkdir -p $folder_genomes
cd "$folder_shell/shell_scripts"
