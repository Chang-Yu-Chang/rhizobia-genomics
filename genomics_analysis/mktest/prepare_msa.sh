#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# This script prepares the MSAs for MK test. It includes the following steps
# 1. Prepare BLAST db from reference genomes, which will be the outgroup. These genomes should have been annotated using Prokka so GFF available
# 2. Move the first sequence from each of the core gene MSA into one single fasta
# 3. Blast the fasta (the number of sequences in this file should be identical to the number of core genes)
# 4. For each gene with a match with reference genome, align the reference sequence to the MSA.
# This reference is the outgroup and should be positioned as the last sequence for the python script to work

# Prepare BLAST database from three genomes
# S. meliloti EM1021
# S. medicae WSM419
# S. fredii NGR234
mkdir -p $folder_data/genomics_analysis/mktest/outgroup/

for ref in em1021 wsm419 ngr234; do

    mkdir -p $folder_data/genomics_analysis/mktest/outgroup/$ref

    # Extract CDS sequences from gff files. Genome fasta are required
    mamba activate gffread
    gffread \
        -x $folder_data/genomics_analysis/mktest/outgroup/$ref/$ref.fasta \
        -g $folder_genomics/fasta/ncbi/$ref.fasta \
        $folder_genomics/gff/$ref.gff
    # -x $folder_data/genomics_analysis/mktest/em1021_genes.fasta denotes the output

    # Make BLAST db
    mamba activate blast
    makeblastdb \
        -in $folder_data/genomics_analysis/mktest/outgroup/$ref/$ref.fasta \
        -dbtype nucl \
        -out $folder_data/genomics_analysis/mktest/outgroup/$ref/$ref
done

# Save the first sequence from each MSA to a fasta
for set_name in elev_med urbn_mel; do
    # #set_name="elev_med"
    # set_name="urbn_mel"

    mkdir -p $folder_data/genomics_analysis/mktest/$set_name/

    # Combined fasta of the first sequence of each core MSA
    #rm $output_fasta
    output_fasta=$folder_data/genomics_analysis/mktest/$set_name/one_seq_from_msa.fasta
    touch $output_fasta
    counter=0

    # Loop over each MSA file in the directory
    for fasta_file in $folder_genomics/pangenome/$set_name/aligned_gene_sequences/*.aln.fas; do
        # Get the gene name from the filename (without the path and extension)
        gene_name=$(basename $fasta_file .aln.fas)

        # Extract the first sequence (header and sequence)
        header_sequence=$(awk '/^>/{header=$0; getline; while(getline && !/^>/) { seq = seq $0} print header; print seq; exit}' $fasta_file)
        sequence=$(echo -e $header_sequence | awk 'NR > 1')

        # Save the extracted sequence under the new gene name
        {
            echo ">$gene_name"  # New header
            echo $sequence  # Sequence without the original header
        } >> $output_fasta

        counter=$(($counter + 1))
        echo $counter
        echo $gene_name
    done

    # Remove the gaps in the sequence that cannot be used in blast
    mamba activate seqkit
    query_fasta=$folder_data/genomics_analysis/mktest/$set_name/one_seq_from_msa_gaps_removed.fasta
    seqkit seq --remove-gaps $output_fasta -o $query_fasta
done


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

# For each gene with a match with reference genome, align the reference sequence to the MSA
function process_msa() {
    local set_name="$1"
    local ref="$2"

    local folder_sourcemsa=$folder_genomics/pangenome/$set_name/aligned_gene_sequences # Source MSA folder
    local folder_outputmsa=$folder_data/genomics_analysis/mktest/$set_name/$ref/msa_with_outgroup # Aligned MSA folder
    local folder_refseq=$folder_data/genomics_analysis/mktest/$set_name/$ref/refseq # Temp folder for extracted sequence from reference
    local blast_results=$folder_data/genomics_analysis/mktest/$set_name/$ref/blast_results.txt # Blast table
    local ref_fasta=$folder_data/genomics_analysis/mktest/outgroup/$ref/$ref.fasta # Reference genome

    mkdir -p $folder_outputmsa
    mkdir -p $folder_refseq

    local counter=0
    while IFS=$'\n' read -r line; do
        # Extract the query name and reference ID using awk
        local query_name=$(echo "$line" | awk -F'\t' '{print $1}')
        local ref_id=$(echo "$line" | awk -F'\t' '{print $2}')

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
    done < $blast_results
}

process_msa "elev_med" "em1021"
process_msa "elev_med" "ngr234"
process_msa "urbn_mel" "wsm419"
process_msa "urbn_mel" "ngr234"














