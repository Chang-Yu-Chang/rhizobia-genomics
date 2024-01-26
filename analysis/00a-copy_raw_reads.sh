#!/usr/bin/env zsh
source ~/.zshrc

folder_genomics=$folder_data/genomics
table_file=$folder_data/temp/00-genomes.csv

# Mapping files for genomics data
table_file=$folder_data/temp/00-genomes.csv
batch_names=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 1)}")
sample_names=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 2)}")
genome_ids=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 3)}")

for i in {1..38}
do 
    cp -u $folder_data/raw/$batch_names[$i]/$sample_names[$i]/reads/raw_reads.fastq.gz \
        $folder_genomics/raw_reads/$genome_ids[$i].fastq.gz
done

