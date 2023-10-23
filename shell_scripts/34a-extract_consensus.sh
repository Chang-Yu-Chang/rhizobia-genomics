cd
source ~/.zshrc

mamba activate bioawk
#mamba env list

# Raw reads
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
folder_raw_result="/Users/cychang/Dropbox/lab/local-adaptation/data/raw/$1"
folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/$1"
sample_id=$2

sample_id=$2

echo $sample_id

consensus_fas="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
consensus_tab="$folder_temp_result/$sample_id/01-filtlong/raw_reads.txt"
# raw_reads="/Users/cychang/Dropbox/lab/local-adaptation/data/raw/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz"
# raw_phred="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/raw_phred.txt"

bioawk -c fastx '{print $name, meanqual($qual), length($seq)}' $raw_reads > $raw_phred
