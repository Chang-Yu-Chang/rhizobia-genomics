#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script assigns taxonomy to the de novo assembled genomes and contigs

cd $folder_shell
refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"
gtdb_db="/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip"
refseq_16s_db="/Users/cychang/bioinformatics/16s/refseq_16s.fasta"

for i in {1..41}
do
    genome_fa="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"
    mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy"

    # Estimate genome distance via mash
    mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/mash"
    zsh 04a-mash.sh \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/mash" \
        $refseq_db

    # Compare genomes via sourmash
    mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/sourmash"
    zsh 04b-sourmash.sh \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/sourmash" \
        $gtdb_db

    # Extract 16S rRNA from genome and blast
    mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s"
    zsh 04c-blast_16s.sh \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/rrna.fasta" \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/rrna.txt" \
        $refseq_16s_db \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/blast.txt" \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/taxonomy.txt"



done

