cd
source ~/.zshrc

# # Install bioawk v1.0
# mamba create -n bioawk
# mamba activate bioawk
# mamba install --yes -c bioconda bioawk=1.0

mamba activate bioawk
#mamba env list

# Raw reads
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
folder_raw_result="/Users/cychang/Dropbox/lab/local-adaptation/data/raw/$1"
folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/"
sample_id=$2

echo $sample_id

raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
raw_phred="$folder_temp_result/$sample_id/01-filtlong/raw_reads.txt"
# raw_reads="/Users/cychang/Dropbox/lab/local-adaptation/data/raw/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz"
# raw_phred="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/raw_phred.txt"

#bioawk -c fastx '{print $qual, length($seq)}' $raw_reads > $raw_phred
bioawk -c fastx '{print $name, $qual, length($seq)}' $raw_reads > $raw_phred

#bioawk -c fastx '{print $name, meanqual($qual), length($seq)} NR==10{exit}' $raw_reads
# # After the first filter
# filtered_reads="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_reads.fastq.gz"
# filtered_phred="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_phred.txt"
#
# bioawk -c fastx '{print $name, meanqual($qual), length($seq)}' $filtered_reads > $filtered_phred
#
# # After the second downsample
# downsampled2_reads="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/03-downsampled2_reads.fastq"
# downsampled2_phred="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/03-downsampled2_phred.txt"
#
# bioawk -c fastx '{print $name, meanqual($qual), length($seq)}' $downsampled2_reads > $downsampled2_phred
