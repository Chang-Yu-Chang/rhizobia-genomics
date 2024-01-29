#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script assigns taxonomy to the de novo assembled genomes and contigs

cd $folder_shell
refseq_db=/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh
gtdb_db=/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip
refseq_16s_db=/Users/cychang/bioinformatics/16s/refseq_16s.fasta
blast_genomes_db=$folder_genomics/blast_db/genomes

# Make customized genome database
zsh 04a-make_db.sh $folder_genomics

mkdir -p $folder_genomics/taxonomy

for i in {1..38}
do
    echo $genome_ids[$i]
    genome_fa=$folder_genomics/genomes/$genome_ids[$i].fasta
    dir=$folder_genomics/taxonomy/$genome_ids[$i]
    mkdir -p $dir

    # # Estimate genome distance via mash
    # mkdir -p $dir/mash
    # zsh 04b-mash.sh \
    #     $genome_fa \
    #     $dir//mash" \
    #     $refseq_db
    
    # Compare strains to database via sourmash
    # mkdir -p $dir/sourmash
    # zsh 04c-sourmash.sh \
    #     $genome_fa \
    #     $dir/sourmash \
    #     $gtdb_db
    
    # Extract 16S rRNA from genome and blast
    mkdir -p $dir/16s
    zsh 04d-blast_16s.sh \
        $genome_fa \
        $dir/16s/rrna.fasta \
        $dir/16s/rrna.txt \
        $refseq_16s_db \
        $dir/16s/blast_16s.txt 

    # Blast genomes to a customized database of meliloti and medicae
    mkdir -p $dir/blast_genome
    zsh 04e-blast_genome.sh \
        $genome_fa \
        $blast_genomes_db \
        $dir/blast_genome/blast_genome.txt

done























