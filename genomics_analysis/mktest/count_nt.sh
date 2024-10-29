#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script counts the nucleotide substitution and generate the contingency tables

mamba activate mktest

function count_nt() {
    local set_name=$1
    local ref=$2

    # Define folder paths
    local folder_msa=$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup_trimmed  # Aligned MSA folder
    local folder_tables=$folder_data/genomics_analysis/mktest/$set_name/$ref/tables
    mkdir -p $folder_tables

    # Process each FASTA file in the folder
    for input_file in "$folder_msa"/*.fasta; do
        local base_name=$(basename "$input_file")
        local gene_name=${base_name%.fasta}

        python sfs_from_fasta_2.py \
            multiFasta \
            $folder_msa/$gene_name.fasta \
            $gene_name.daf \
            $gene_name.div \
            standard \
            "$folder_tables/"

        echo "Finished: $gene_name"
    done
}

count_nt "elev_med" "em1021"
count_nt "elev_med" "ngr234"
count_nt "urbn_mel" "wsm419"
count_nt "urbn_mel" "ngr234"

