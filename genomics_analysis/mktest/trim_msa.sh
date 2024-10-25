#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh


function trim_msa() {
    local set_name="$1"
    local ref="$2"

    # Define input and output directory paths
    local input_directory="$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup"  # Aligned MSA folder
    local output_directory="$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup_trimmed"  # Trimmed

    # Create the output directory if it doesn't exist
    mkdir -p "$output_directory"

    mamba activate seqkit

    # Process each FASTA file in the input directory
    for input_file in "$input_directory"/*.fasta; do
        # Get the base name for the output file
        local base_name=$(basename "$input_file")
        local output_file="$output_directory/$base_name"

        echo "Processing $base_name..."

        # Use seqkit to check if sequences start with "ATG"
        if seqkit grep -s -p ATG "$input_file" | seqkit stat | grep -q '0 sequences'; then
            echo "Discarding $base_name: not all sequences start with 'ATG'."
        else
            # Trim sequences to the nearest length divisible by 3
            awk -v output="$output_file" '
                /^>/ { print > output; header = $0; next }
                {
                    seq_length = length($0);
                    trimmed_length = seq_length - (seq_length % 3);  # Trim to nearest divisible by 3
                    print substr($0, 1, trimmed_length) > output;     # Print the trimmed sequence
                }
            ' "$input_file"

            echo "Saved processed sequences to: $output_file"
        fi

    done
}

#trim_msa "elev_med" "em1021"
trim_msa "elev_med" "ngr234"
#trim_msa "urbn_mel" "wsm419"
trim_msa "urbn_mel" "ngr234"




# for set_name in elev_med urbn_mel; do
#     input_directory=$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup  # aligned MSA folder
#     output_directory=$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup_trimmed  # trimmed
#
#     # Create the output directory if it doesn't exist
#     mkdir -p "$output_directory"
#
#     mamba activate seqkit
#     #base_name=AAH1.fasta
#     # Process each FASTA file in the input directory
#     for input_file in $input_directory/*.fasta; do
#         # Get the base name for the output file
#         base_name=$(basename $input_file)
#         output_file=$output_directory/$base_name
#
#         echo "Processing $base_name..."
#
#         # Use seqkit to get sequences that do not start with "ATG"
#         if seqkit grep -s -p ATG $input_file | seqkit stat | grep -q '0 sequences'; then
#             # If no sequences start with "ATG", skip this file
#             echo "Discarding $base_name: not all sequences start with 'ATG'."
#         else
#             # Trim sequences to the nearest length divisible by 3
#             awk -v output="$output_file" '
#                 /^>/ { print > output; header = $0; next }
#                 {
#                     seq_length = length($0);
#                     trimmed_length = seq_length - (seq_length % 3);  # Trim to nearest divisible by 3
#                     print substr($0, 1, trimmed_length) > output;     # Print the trimmed sequence
#                 }
#             ' "$input_file"
#
#             echo "Saved processed sequences to: $output_file"
#         fi
#
#     done
#
# done
