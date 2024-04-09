#!/usr/bin/env zsh
source ~/.zshrc

# This blasts the whole genome across the database
# $1: genome in fasta
# $2: genomes_db
# $3: output blast prediction in txt

# Blast
mamba activate blast

# Concatenate genome
blastn -query $1 -db $2 -out $3 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 5
#blastn -query $1 -db $2 -out $3 -outfmt '6 qseqid sseqid stitle bitscore evalue length pident' -num_alignments 10
#blastn -query your_query_sequence.fasta -db "$folder_genomics/blast_db/genomes" -out "$folder_genomics/blast_db/genomes"
