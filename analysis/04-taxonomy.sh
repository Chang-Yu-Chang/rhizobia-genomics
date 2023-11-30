#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script assigns taxonomy to the de novo assembled genomes and contigs

cd $folder_shell
echo "04-taxonomy"
refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"
gtdb_db="/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip"

for i in {2..19}
do
    echo "$folder_raw/$batch_ids[$i]/$sample_ids[$i]"
    genome_fa="$folder_genomes/$sample_ids[$i]/03-denovo_assembly/genome.fasta"

    # Estimate genome and metagenome distance via mash
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



done

