#!/usr/bin/env zsh
source ~/.zshrc

# This prepares the anvio object `contigs_db` for each genome
# $1 parser directory
# $2: prokka annotation in gff
# $3: prokka gene calls in txt
# $4: prokka gene annot in txt
# $5: target genome in fasta
# $6: target genome in anvio contig db
# $7: sample id
# $8: target genome contigs_stat in txt

conda activate
mamba activate anvio-8

# Import prokka annotations into anvio. Prokka should have had happened
python3 $1 $2 --gene-calls $3 --annotation $4

# Create the contig database attached with external gene call
anvi-gen-contigs-database \
    --force-overwrite \
    -f $5 \
    -o $6 \
    --external-gene-calls $3 \
    --project-name $7 \
    --ignore-internal-stop-codons
# --external-gene-calls A TAB-delimited file to define external gene calls.
# --ignore-internal-stop-codons This is only relevant when you have an external gene calls file

# Import the functional annotations
anvi-import-functions \
    -c $6 \
    -i $4

# Run hmms, which helps annotate the genes in your contigs-db. This adds the HMM info to the database
anvi-run-hmms -c $6

# Output the contig stat from each contigs-db
anvi-display-contigs-stats \
    $6 \
    --report-as-text \
    -o $8
# --report-as-text  If you give this flag, Anvi'o will not open new browser to show Contigs database statistics and write all stats to TAB separated file and you should also give --output-file with this flag otherwise Anvi'o will complain.












