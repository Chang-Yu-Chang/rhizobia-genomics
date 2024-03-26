#!/usr/bin/env zsh
source ~/.zshrc
source ../analysis/00-env_vars.sh


# for i in Chang_Q5C_{1..10} Chang_Q5C_{12..17} Chang_Q5C_19
# do
#     cp "$folder_raw/Chang_Q5C_results/$i/reads/raw_reads.fastq.gz" \
#         "$folder_genomics/raw_reads/$i.fastq.gz"
# done

for i in Chang_Q5C_11 Chang_Q5C_18
do
    cp "$folder_raw/Chang_Q5C_results_repeated/$i/reads/raw_reads.fastq.gz" \
        "$folder_genomics/raw_reads/$i.fastq.gz"
done

for i in Chang_W8S_{1..18}
do
    cp "$folder_raw/Chang_W8S_results/$i/reads/raw_reads.fastq.gz" \
        "$folder_genomics/raw_reads/$i.fastq.gz"
done
