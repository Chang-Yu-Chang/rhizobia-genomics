#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assigns taxonomy to genomes using sourmash

cd $folder_shell
gtdb_db=/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip

mkdir -p $folder_genomics/taxonomy

# Compare strains to database via sourmash

for i in {30..32}
do
    echo $genome_ids[$i]
    genome_fa=$folder_genomics/fasta/genomes/$genome_ids[$i].fasta
    dir=$folder_genomics/taxonomy/$genome_ids[$i]

    mkdir -p $dir/sourmash
    zsh 04c-sourmash.sh \
        $genome_fa \
        $dir/sourmash \
        $gtdb_db
done

