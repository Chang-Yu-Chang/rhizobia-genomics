import os
import argparse
from Bio import SeqIO

def read_genomes_list(genomes_file):
    with open(genomes_file, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def get_gene_name(filename):
    # Derive gene name from filename (modify if needed)
    return os.path.splitext(os.path.basename(filename))[0]

def count_snps(sequences):
    seqs = [str(record.seq) for record in sequences]
    seq_length = len(seqs[0])
    snp_count = 0
    for i in range(seq_length):
        column = [seq[i] for seq in seqs]
        if len(set(column)) > 1:
            snp_count += 1
    return snp_count, seq_length

def process_folder(folder_path, target_genomes):
    gene_results = []  # List of (gene_name, gene_length, snp_count)
    for filename in os.listdir(folder_path):
        if filename.endswith('.fas') or filename.endswith('.fasta'):
            filepath = os.path.join(folder_path, filename)
            sequences = list(SeqIO.parse(filepath, 'fasta'))

            # Filter sequences to target genomes
            filtered_sequences = []
            for record in sequences:
                header = record.id
                genome_id = header.split(';')[0]
                if genome_id in target_genomes:
                    filtered_sequences.append(record)
            
            # Check if all target genomes are present exactly once
            present_genomes = set(r.id.split(';')[0] for r in filtered_sequences)
            if present_genomes == target_genomes and len(filtered_sequences) == len(target_genomes):
                snp_num, gene_length = count_snps(filtered_sequences)
                gene_name = get_gene_name(filename)
                gene_results.append((gene_name, gene_length, snp_num))
            else:
                print(f"Skipping {filename}: not all target genomes present exactly once.")
    return gene_results

def main():
    parser = argparse.ArgumentParser(description='Count SNPs and report gene lengths among a subset of genomes.')
    parser.add_argument('folder', type=str, help='Path to folder containing alignment files.')
    parser.add_argument('genomes_list', type=str, help='Path to file with target genome IDs (one per line).')
    parser.add_argument('output', type=str, help='Path to output txt file with results.')
    args = parser.parse_args()

    target_genomes = read_genomes_list(args.genomes_list)

    results = process_folder(args.folder, target_genomes)

    # Save results to output file
    with open(args.output, 'w') as out_f:
        out_f.write("GeneName\tGeneLength\tNumSNPs\n")
        for gene_name, gene_length, snp_count in results:
            out_f.write(f"{gene_name}\t{gene_length}\t{snp_count}\n")
    
    print(f"Results saved to {args.output}")

if __name__ == '__main__':
    main()
