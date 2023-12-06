#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script prepares the dataset for pangenome analysis

cd $folder_shell
echo "06-pangenome_prep"
gff_parser=~/bioinformatics/anvio/gff_parser.py # For using the gff_parser.py. This has to be downloaded

for i in {24..41}
do
    genome_fa="$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"

    # Prepare anvio dataset
    zsh 06a-prepare_anvio.sh \
        $gff_parser \
        "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka/annotated.gff" \
        "$folder_genomes/$sample_ids[$i]/06-pangenome_prep/gene_calls.txt" \
        "$folder_genomes/$sample_ids[$i]/06-pangenome_prep/gene_annot.txt" \
        $genome_fa \
        "$folder_genomes/$sample_ids[$i]/06-pangenome_prep/genome.db" \
        $sample_ids[$i] \
        "$folder_genomes/$sample_ids[$i]/06-pangenome_prep/genome_stat.txt"

    # Annotate genomes via prokka
    # mkdir -p "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka"
    # zsh 06b-anvio_prep.sh \
    #     $genome_fa \
    #     "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka"

done
