cd
source ~/.zshrc

# Polish Flye assembly bia Medaka
mamba activate medaka
mamba env list

raw_reads=$1
draft_assembly=$2
medaka_folder=$3


# raw_reads=~/Dropbox/lab/local-adaptation/data/raw/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz
# draft_assembly=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/assembly.fasta
# medaka_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka

medaka_consensus -i $raw_reads -d $draft_assembly -o $medaka_folder -t 10 -m r941_min_high_g303
# `-i "$raw_reads"` fastx input basecalls (required).
# `-d "$draft_assembly"` fasta input assembly (required).
# `-o "$medaka_folder"` output folder (default: medaka).
# `-t 10` number of threads with which to create features (default: 1).
# `-m r941_min_high_g303` medaka model, (default: r941_min_high_g360).
