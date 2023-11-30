#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script annotates the de novo assembled genomes

cd $folder_shell
echo "05-gene_annotation"
bakta_db="/Users/cychang/bioinformatics/bakta/db"

for i in {2..19}
do
    #echo "$folder_raw/$batch_names[$i]/$sample_ids[$i]"
    genome_fa="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"

    # Annotate genomes via prokka
    mkdir -p "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka"
    zsh 05a-prokka.sh \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka"

    # Annotate genomes via bakta
    mkdir -p "$folder_genomes/$sample_ids[$i]/05-gene_annotation/bakta"
    zsh 05b-bakta.sh \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/05-gene_annotation/bakta" \
        $bakta_db # This database is mandatory and must be downloaded before annotation

done
