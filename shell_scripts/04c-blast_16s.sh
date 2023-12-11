#!/usr/bin/env zsh
source ~/.zshrc

# This extracts the rRNA sequences from the genome fasta and blast
# $1: genome in fasta
# $2: output rrna in fasta
# $3: output rrna in a table in txt
# $4: refseq_16s_db
# $5: output blast prediction in txt
# $6: output taxonomy prediction in txt

# Extracrt rRNA
conda activate
mamba activate barrnap

barrnap -o $2 < $1 > $3


# Blast
conda activate
mamba activate blast

blastn -query $2 -db $4 -out $5 -outfmt '6 qseqid sseqid stitle bitscore evalue length pident' -num_alignments 10
