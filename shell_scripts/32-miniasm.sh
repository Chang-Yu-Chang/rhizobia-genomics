cd
source ~/.zshrc


# Downsample the reads to 250 Mbp
mamba activate filtlong
mamba env list

filtered_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_reads.fastq.gz
downsample_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/02-downsampled_reads.fastq

filtlong --target_bases 250000000 "$filtered_reads" > "$downsample_reads"
# `--target_bases 250000000` remove the worst reads unitl only 250 Mbp remain

# Assembly
mamba activate miniasm
mamba env list

assemble_paf=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/02-miniasm_pilot.paf.gz
assemble_gfa=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/02-miniasm_pilot.gfa
assemble_fq=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/02-miniasm_pilot.fasta

## Overlap raw reads using minimap
minimap2 -x ava-ont -t8 "$downsample_reads" "$downsample_reads" | gzip -1 > "$assemble_paf"
# `-x ava-ont` present Nanopore read overlap
# `-t8` use 8 threads
# `gzip -1` fastest compression

## Assemble reads using miniams
miniasm -f "$downsample_reads" "$assemble_paf" > "$assemble_gfa"
# `-f` read sequences

## Convert assembly graph to fastq
any2fasta "$assemble_gfa" > "$assemble_fq"

