#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script prepares the dataset for pangenome analysis

cd $folder_shell
echo "06-pangenome_prep"
#bakta_db="/Users/cychang/bioinformatics/bakta/db"

for i in {1..19}
do
    echo "$folder_raw/$batch_ids[$i]/$sample_ids[$i]"
    genome_fa="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"

    # Prepare anvio dataset



    # Annotate genomes via prokka
    mkdir -p "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka"
    zsh 06b-anvio_prep.sh \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka"

done
