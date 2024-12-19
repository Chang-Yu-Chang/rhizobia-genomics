#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assigns taxonomy to the de novo assembled genomes and contigs

cd $folder_shell
refseq_16s_db=/Users/cychang/bioinformatics/16s/refseq_16s.fasta
blast_genomes_db=$folder_genomics/blast_db/genomes

mkdir -p $folder_genomics/taxonomy

for i in {1..36}
do
    echo $genome_ids[$i]
    genome_fa=$folder_genomics/fasta/genomes/$genome_ids[$i].fasta
    dir=$folder_genomics/taxonomy/$genome_ids[$i]
    mkdir -p $dir

    mkdir -p $dir/16s
    zsh 04d-blast_16s.sh \
        $genome_fa \
        $dir/16s/rrna.fasta \
        $dir/16s/rrna.txt \
        $refseq_16s_db \
        $dir/16s/blast_16s.txt

    # Blast genomes to a customized database of ensifer strains
    mkdir -p $dir/blast_genome
    zsh 04e-blast_genome.sh \
        $genome_fa \
        $blast_genomes_db \
        $dir/blast_genome/blast_genome.txt

done
