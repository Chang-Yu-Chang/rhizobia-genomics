#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script prepares the MSAs for MK test. It includes the following steps

mamba activate mktest

set_name="elev_med"
# set_name="urbn_mel"

for set_name in elev_med urbn_mel; do

    folder_msa=$folder_data/genomics_analysis/mktest/$set_name/msa_with_outgroup_trimmed  # aligned MSA folder
    folder_tables=$folder_data/genomics_analysis/mktest/$set_name/tables
    mkdir -p $folder_tables

    #gene_name=accA
    for input_file in $folder_msa/*.fasta; do
        base_name=$(basename $input_file)
        gene_name=${base_name%.fasta}

        python sfs_from_fasta_2.py \
            multiFasta \
            $folder_msa/$gene_name.fasta \
        	$gene_name.daf \
        	$gene_name.div \
        	standard \
        	$folder_tables/

        echo "Finished: $gene_name"
    done
done
