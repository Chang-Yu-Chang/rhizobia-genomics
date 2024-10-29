#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script prepares the MSAs for MK test. It includes the following steps
# 1. Prepare BLAST db from reference genomes, which will be the outgroup. These genomes should have been annotated using Prokka so GFF available
# 2. Move the first sequence from each of the core gene MSA into one single fasta `one_seq_from_msa`. Also remove the gaps since it's not allowed in blast
# 3. Blast the query fasta sequences (the number of sequences in this file should be identical to the number of core genes)
# 4. Distinct the blast results. Some gene has multiple blast hits. Use only the top hit per query gene
# 5. For each gene with a match with reference genome, align the reference sequence to the MSA. The sequences from the outgroup/reference genome are stored in ref_seq/
# 6. Trims the MSA to be divisible by 3, and filter for those CDS starting with ATG
# This reference is the outgroup and should be positioned as the last sequence for the python script to work

# Prepare BLAST database from three genomes
# S. meliloti EM1021
# S. medicae WSM419
# S. fredii NGR234
mkdir -p $folder_data/genomics_analysis/mktest/outgroup/
make_blast_db() {
    local ref="$1"

    # Create the necessary directory
    mkdir -p $folder_data/genomics_analysis/mktest/outgroup/$ref

    # Extract CDS sequences from GFF files. Genome FASTA are required
    mamba activate gffread
    gffread \
        -x $folder_data/genomics_analysis/mktest/outgroup/$ref/$ref.fasta\
        -g $folder_genomics/fasta/ncbi/$ref.fasta\
        $folder_genomics/gff/$ref.gff

    # Make BLAST db
    mamba activate blast
    makeblastdb \
        -in $folder_data/genomics_analysis/mktest/outgroup/$ref/$ref.fasta\
        -dbtype nucl \
        -out $folder_data/genomics_analysis/mktest/outgroup/$ref/$ref
}

make_blast_db em1021
make_blast_db wsm419
make_blast_db ngr234

# Save the first sequence from each MSA to a fasta
save_first_seq() {
    local set_name="$1"
    mkdir -p $folder_data/genomics_analysis/mktest/$set_name/

    # Combined fasta of the first sequence of each core MSA
    local output_fasta=$folder_data/genomics_analysis/mktest/$set_name/one_seq_from_msa.fasta

    # Remove the output_fasta if it exists
    if [[ -f $output_fasta ]]; then
        rm $output_fasta
    fi

    touch $output_fasta
    local counter=0

    # Loop over each MSA file in the directory
    for fasta_file in $folder_genomics/pangenome/$set_name/aligned_gene_sequences/*.aln.fas; do
        # Get the gene name from the filename (without the path and extension)
        local gene_name=$(basename $fasta_file .aln.fas)

        # Extract the first sequence (header and sequence)
        local header_sequence=$(awk '/^>/{header=$0; getline; while(getline && !/^>/) { seq = seq $0} print header; print seq; exit}' "$fasta_file")
        local sequence=$(echo -e "$header_sequence" | awk 'NR > 1')

        # Save the extracted sequence under the new gene name
        {
            echo ">$gene_name"  # New header
            echo $sequence   # Sequence without the original header
        } >> $output_fasta

        counter=$((counter + 1))
        echo $counter
        echo $gene_name
    done

    # Remove the gaps in the sequence that cannot be used in blast
    mamba activate seqkit
    local query_fasta=$folder_data/genomics_analysis/mktest/$set_name/one_seq_from_msa_gaps_removed.fasta
    seqkit seq --remove-gaps $output_fasta -o $query_fasta
}

save_first_seq "elev_med"
save_first_seq "urbn_mel"

# Blast the query against the reference
function run_blast() {
    local set_name="$1"
    local ref="$2"
    local query_fasta=$folder_data/genomics_analysis/mktest/$set_name/one_seq_from_msa_gaps_removed.fasta
    local blast_db=$folder_data/genomics_analysis/mktest/outgroup/$ref/$ref

    mkdir -p $folder_data/genomics_analysis/mktest/$set_name/$ref/

    mamba activate blast
    blastn \
        -query $query_fasta\
        -db $blast_db\
        -out $folder_data/genomics_analysis/mktest/$set_name/$ref/blast_results.txt\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -num_alignments 5
}

run_blast "elev_med" "em1021"
run_blast "elev_med" "ngr234"
run_blast "urbn_mel" "wsm419"
run_blast "urbn_mel" "ngr234"

# Some gene has multiple blast hits. Use only the top hit per query gene
function distinct_blast() {
    local set_name="$1"
    local ref="$2"

    Rscript distinct_blast.R $set_name $ref \
        $folder_data/genomics_analysis/mktest/$set_name/$ref/blast_results.txt \
        $folder_data/genomics_analysis/mktest/$set_name/$ref/blast_distinct.txt
}

distinct_blast "elev_med" "em1021"
distinct_blast "elev_med" "ngr234"
distinct_blast "urbn_mel" "wsm419"
distinct_blast "urbn_mel" "ngr234"

# For each gene with a match with reference genome, align the reference sequence to the MSA
function process_msa() {
    local set_name="$1"
    local ref="$2"

    local list_sccg=$folder_data/genomics_analysis/gene_content/$set_name/list_sccg.csv # The list of single copy core genes
    local folder_sourcemsa=$folder_genomics/pangenome/$set_name/aligned_gene_sequences # Source MSA folder
    local folder_outputmsa=$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup # Aligned MSA folder
    local folder_refseq=$folder_data/genomics_analysis/mktest/$set_name/$ref/refseq # Temp folder for extracted sequence from reference
    local blast_results=$folder_data/genomics_analysis/mktest/$set_name/$ref/blast_results.txt # Blast table
    local ref_fasta=$folder_data/genomics_analysis/mktest/outgroup/$ref/$ref.fasta # Reference genome

    mkdir -p $folder_outputmsa
    mkdir -p $folder_refseq

    # Read single-copy core genes into an array
    local single_copy_genes=()
    while IFS=',' read -r gene; do
        single_copy_genes+=("$gene")  # Add gene to array
    done < $list_sccg

    echo $single_copy_genes

    local counter=0
    while IFS=$'\n' read -r line; do
        # Extract the query name and reference ID using awk
        local query_name=$(echo "$line" | awk -F'\t' '{print $1}')
        local ref_id=$(echo "$line" | awk -F'\t' '{print $2}')

        # Check if query_name exists in the single_copy_genes array
        if [[ " ${single_copy_genes[*]} " == *" $query_name "* ]]; then
            # Extract the matched reference sequence using seqkit
            mamba activate seqkit
            local ref_seq=$(seqkit grep -n -p $ref_id $ref_fasta)  # Adjust based on how ref_id is formatted

            # Create a temporary FASTA file for the reference sequence
            echo $ref_seq> $folder_refseq/$query_name.fasta

            # Find corresponding MSA file in the specified directory
            local msa_file=$folder_sourcemsa/${query_name}.aln.fas # Assuming MSAs are named after query names

            # Use MAFFT to align the matched reference sequence to the corresponding MSA
            mamba activate mafft
            mafft --add $folder_refseq/$query_name.fasta $msa_file> $folder_outputmsa/${query_name}.fasta

            counter=$((counter + 1))
            echo $counter
            echo $query_name
        else
            echo "Skipping $query_name: not a single-copy core gene."
        fi
    done < $blast_results
}

process_msa "elev_med" "em1021"
process_msa "elev_med" "ngr234"
process_msa "urbn_mel" "wsm419"
process_msa "urbn_mel" "ngr234"

# Trims the MSA to be divisible by 3, and filter for those CDS starting with ATG
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

trim_msa "elev_med" "em1021"
trim_msa "elev_med" "ngr234"
trim_msa "urbn_mel" "wsm419"
trim_msa "urbn_mel" "ngr234"
