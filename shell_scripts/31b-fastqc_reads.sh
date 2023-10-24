cd
source ~/.zshrc

mamba activate fastqc
#mamba activate multiqc

raw_reads=$1
filtered_reads=$2

raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
filtered_reads="$folder_temp_result/$sample_id/01-filtlong/01-filtered_reads.fastq.gz"

raw_reads="$folder_raw_result/$sample_id/reads/raw_reads.fastq.gz"
filtered_reads="$folder_temp_result/$sample_id/01-filtlong/01-filtered_reads.fastq.gz"


fastqc $filtered_reads
