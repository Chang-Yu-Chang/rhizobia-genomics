#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script assigns taxonomy the de novo assembled genomes

cd $folder_shell
echo "04-taxonomy"

for i in {1..19}
do
    echo "$folder_raw/$batch_ids[$i]/$sample_ids[$i]"
    medaka_consensus="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/medaka/consensus.fasta"
    # raw_reads="$folder_raw/$batch_ids[$i]/$sample_ids[$i]/reads/raw_reads.fastq.gz"
    # filtered_reads="$folder_genomes/$sample_ids[$i]/01-reads_qc/filtered_reads.fastq.gz"
    # downsampled_reads="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/downsampled_reads.fastq"
    # assembly_prefix="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/miniasm/draft_assembly"
    # draft_genome="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/miniasm/draft_genome.fasta"
    # downsampled_reads2="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/downsampled_reads2.fastq"
    refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"
    gtdb_db="/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip"

    # Estimate genome and metagenome distance via mash
    mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/mash"
    zsh 04a-mash.sh \
        $medaka_consensus \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/mash" \
        $refseq_db

    # Compare genomes via sourmash
    mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/sourmash"
    zsh 04b-sourmash.sh \
        $medaka_consensus \
        "$folder_genomes/$sample_ids[$i]/04-taxonomy/sourmash" \
        $gtdb_db



done

