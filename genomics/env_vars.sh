#!/usr/bin/env zsh
source ~/.zshrc

# This script stores the environment varaibles shared by all shell scripts

folder_shell=/Users/cychang/Desktop/lab/rhizobia-genomics/genomics/genomics_cli
folder_data=/Users/cychang/Dropbox/lab/rhizobia-genomics/data
folder_raw=$folder_data/raw
folder_temp=$folder_data/temp
folder_genomics=$folder_data/genomics

# Mapping files for genomics data
table_file=$folder_data/mapping/genomes.csv
batch_names=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 1)}")
genome_names=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 2)}")
genome_ids=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 3)}")
