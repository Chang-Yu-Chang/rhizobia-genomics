import os
import pandas as pd
import argparse
from Bio import SeqIO

def concatenate_genomes(input_directory, output_file, genes_csv_file, genomes_csv_file):
    # Read the gene names from the CSV file
    genes_df = pd.read_csv(genes_csv_file, header=None)
    gene_names = genes_df[0].tolist()  # List of gene names

    # Read the genome IDs from the CSV file
    genomes_df = pd.read_csv(genomes_csv_file, header=None)
    genome_ids = genomes_df[0].tolist()  # List of genome IDs

    # Initialize a dictionary to hold sequences for each genome ID
    sequences = {genome_id: '' for genome_id in genome_ids}  # Empty string for collecting sequences

    # Loop through all alignment files in the directory (each representing a gene)
    for filename in os.listdir(input_directory):
        if filename.endswith('.aln.fas'):
            gene_name = filename[:-8]  # Extracts gene name from filename (e.g., alaA from alaA.aln.fas)
            if gene_name not in gene_names:
                continue  # Skip this file if it's not in the gene list

            file_path = os.path.join(input_directory, filename)
            
            # Read sequence records from the file
            for record in SeqIO.parse(file_path, 'fasta'):
                 # Extract the genome ID from the header, accounting for leading underscores
                header_parts = record.id.split(';')[0].strip()  # Get the part before the semicolon
                # Check for leading underscore
                genome_id = header_parts.split('_')[-1]  # This will get 'g8' from '_R_g8' or 'g8' from 'g8; someinfo'

                # If genome ID is in our list, concatenate the sequence
                if genome_id in sequences:
                    sequences[genome_id] += str(record.seq)  # Concatenate the sequence

    # Write the concatenated sequences to the output file
    with open(output_file, 'w') as out_file:
        for genome_id, seq in sequences.items():
            # Only write the header and sequence if the sequence is not empty
            if seq:
                out_file.write(f">{genome_id}\n{seq}\n")  # Write the genome ID as the header

    print(f"Combined alignment written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Concatenate genome sequences from multiple gene alignment files.')
    parser.add_argument('input_directory', type=str, help='Path to the folder containing alignment files')
    parser.add_argument('output_file', type=str, help='Name of the output alignment file')
    parser.add_argument('genes_csv_file', type=str, help='Path to the CSV file containing gene names')
    parser.add_argument('genomes_csv_file', type=str, help='Path to the CSV file containing genome IDs')
    
    args = parser.parse_args()

    # Call the function to concatenate genomes
    concatenate_genomes(args.input_directory, args.output_file, args.genes_csv_file, args.genomes_csv_file)

if __name__ == "__main__":
    main()
