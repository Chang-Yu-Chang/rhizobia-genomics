#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script prepares the MSAs for MK test. It includes the following steps

mamba activate mktest

#set_name="elev_med"
# set_name="urbn_mel"

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

count_nt "elev_med" "ngr234"
count_nt "urbn_mel" "ngr234"



# for set_name in elev_med urbn_mel; do
#
#     folder_msa=$folder_data/genomics_analysis/mktest/$set_name/msa_with_outgroup_trimmed  # aligned MSA folder
#     folder_tables=$folder_data/genomics_analysis/mktest/$set_name/tables
#     mkdir -p $folder_tables
#
#     #gene_name=accA
#     for input_file in $folder_msa/*.fasta; do
#         base_name=$(basename $input_file)
#         gene_name=${base_name%.fasta}
#
#         python sfs_from_fasta_2.py \
#             multiFasta \
#             $folder_msa/$gene_name.fasta \
#         	$gene_name.daf \
#         	$gene_name.div \
#         	standard \
#         	$folder_tables/
#
#         echo "Finished: $gene_name"
#     done
# done
