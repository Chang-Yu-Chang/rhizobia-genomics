cd
source ~/.zshrc


cd ~/Deskd
zsh 31

raw_reads=~/Dropbox/lab/local-adaptation/data/raw/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz
filtered_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_reads.fastq.gz
#log_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/logs/01-filter_reads.log

filtlong --keep_percent 95 "$raw_reads" | gzip > "$filtered_reads"
# `--keep_percent 95` throw out the worst 5% of reads
