import os
import pandas as pd
from Bio import SeqIO

# Define the input directory and output file
input_directory = '/Users/cychang/Dropbox/lab/rhizobia-genomics/data/genomics/pangenome/elev_med/aligned_gene_sequences' 
output_file = '/Users/cychang/Dropbox/lab/rhizobia-genomics/data/genomics/pangenome/elev_med/combined_sccg_alignment.fas'
genes_csv_file = '/Users/cychang/Dropbox/lab/rhizobia-genomics/data/genomics_analysis/gene_content/elev_med/list_sccg.csv'

genomes_csv_file = 'path/to/your/genomes.csv'  # Change to your actual genomes CSV file path

# Read the gene names from the CSV file
genes_df = pd.read_csv(genes_csv_file, header=None)
gene_names = genes_df[0].tolist()  # List of gene names

# Read the genome IDs from the CSV file
# genomes_df = pd.read_csv(genomes_csv_file, header=None)
# genome_ids = genomes_df[0].tolist()  # List of genome IDs
genome_ids = ['g4', 'g5', 'g6', 'g8', 'g9', 'g11', 'g13', 'g16', 'g17', 'g19']

# Initialize a dictionary to hold sequences for each genome ID
sequences = {genome_id: '' for genome_id in genome_ids}  # Empty string for collecting sequences

# Loop through all alignment files in the directory (each representing a gene)
for filename in os.listdir(input_directory):
    if filename.endswith('.aln.fas'):
        gene_name = filename[:-8]  # Extracts gene name from filename (e.g., from `alaA.aln.fas`)
        if gene_name not in gene_names:
            continue  # Skip this file if it's not in the gene list

        file_path = os.path.join(input_directory, filename)
        
        # Read sequence records from the file
        for record in SeqIO.parse(file_path, 'fasta'):
            # Assume the record ID has the format 'g1_some_additional_info'
            header_parts = record.id.split(';')[0].strip()  # Get the header part before the semicolon
            genome_id = header_parts.split('_')[0]  # Extract the genome ID

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
