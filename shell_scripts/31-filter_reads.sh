cd
source ~/.zshrc

# Throw out the worst 5% reads
conda activate
mamba activate filtlong
mamba env list

raw_reads=$1
filtered_reads=$2

filtlong --keep_percent 95 $raw_reads | gzip > $filtered_reads
# `--keep_percent 95` throw out the worst 5% of reads
