cd
source ~/.zshrc


# Using information acquired from the Miniasm assembly, re-downsample the reads to ~100x coverage (do nothing if there isn't at least 100x coverage) with heavy weight applied to removing low quality reads (helps small plasmids stick around)
mamba activate filtlong
mamba env list

filtered_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_reads.fastq.gz
assemble_fq=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/02-miniasm_pilot.fasta
downsampled_reads2=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka/03-downsampled2_reads.fastq
genome_size=$(grep -v ">" "$assemble_fq" | grep -E -o "G|C|T|A|N" | wc -l)

filtlong --target_bases $((genome_size * 100)) --mean_q_weight 10 --assembly "$assemble_fq" "$filtered_reads" | > "$downsampled_reads2"
# `--target_bases $((genome_size * 100))` allow the coverage X100 of genome size
# `--mean_q_weight 10` specify a mean quality weight of 10 (instead of the default of 1) makes mean read quality more important when choosing the best reads
# `--assembly "$assemble_fq"` use reference assembly in FASTA format



# Select for high quality ONT reads and assemble again
mamba activate flye
mamba env list

downsampled_reads2=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka/03-downsampled2_reads.fastq
assembly_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka/

flye --meta --nano-corr "$downsampled_reads2" --out-dir "$assembly_folder" --threads 20
# `--meta` metagenome / uneven coverage mode
# `--nano-corr "$downsampled_reads2"` ONT reads that were corrected with other methods (<3% error)
# `--out-dir "$assembly_folder"` Output directory
# `--threads 20`  number of parallel threads [1]


# Polish Flye assembly bia Medaka
mamba activate medaka
mamba env list

raw_reads=~/Dropbox/lab/local-adaptation/data/raw/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz
draft_assembly=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka/assembly.fasta
medaka_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka/medaka/

medaka_consensus -i "$raw_reads" -d "$draft_assembly" -o "$medaka_folder" -t 10 -m r941_min_high_g303
# `-i "$raw_reads"` fastx input basecalls (required).
# `-d "$draft_assembly"` fasta input assembly (required).
# `-o "$medaka_folder"` output folder (default: medaka).
# `-t 10` number of threads with which to create features (default: 1).
# `-m r941_min_high_g303` medaka model, (default: r941_min_high_g360).






















