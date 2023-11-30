#!/usr/bin/env zshs
source ~/.zshrc

# This script saves the environment varaibles shared by all shell scripts

folder_shell="/Users/cychang/Desktop/lab/local-adaptation/shell_scripts"
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
folder_raw="$folder_data/raw"
folder_temp="$folder_data/temp"
folder_genomics="$folder_data/genomics"
folder_genomes="$folder_data/genomics/genomes"

# Mapping files for genomics data
table_file="$folder_data/temp/00-genomes_mapping.csv"
batch_ids=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 1)}")
sample_ids=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 2)}")

# Make folders for temporary genomic files
for i in {1..19}
do
    for wf in 01-reads_qc 02-guided_assembly 03-denovo_assembly 04-taxonomy 05-gene_annotation 06-pangenome_prep
    do
        mkdir -p "$folder_genomes/$sample_ids[i]/$wf"
    done
done

# Set up conda and environments. Uncomment this chunk if not set up
# zsh $folder_shell/00a-setup_conda.sh
# zsh $folder_shell/00b-setup_envs.sh
