cd
source ~/.zshrc

filtered_reads=$1
assemble_fq=$2
downsampled_reads2=$3
genome_size=$4
assembly_folder=$5


# Using information acquired from the Miniasm assembly, re-downsample the reads to ~100x coverage (do nothing if there isn't at least 100x coverage) with heavy weight applied to removing low quality reads (helps small plasmids stick around)
mamba activate filtlong
mamba env list

# filtered_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/01-filtered_reads.fastq.gz
# assemble_fq=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/02-miniasm_pilot.fasta
# downsampled_reads2=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/03-downsampled2_reads.fastq
# genome_size=$(grep -v ">" "$assemble_fq" | grep -E -o "G|C|T|A|N" | wc -l)

filtlong --target_bases $((genome_size * 100)) --mean_q_weight 10 --assembly $assemble_fq $filtered_reads | > $downsampled_reads2
# `--target_bases $((genome_size * 100))` allow the coverage X100 of genome size
# `--mean_q_weight 10` specify a mean quality weight of 10 (instead of the default of 1) makes mean read quality more important when choosing the best reads
# `--assembly "$assemble_fq"` use reference assembly in FASTA format


# Select for high quality ONT reads and assemble again
mamba activate flye
mamba env list

# downsampled_reads2=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/03-downsampled2_reads.fastq
# assembly_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye/

flye --meta --nano-corr $downsampled_reads2 --out-dir $assembly_folder --threads 20
# `--meta` metagenome / uneven coverage mode
# `--nano-corr "$downsampled_reads2"` ONT reads that were corrected with other methods (<3% error)
# `--out-dir "$assembly_folder"` Output directory
# `--threads 20`  number of parallel threads [1]


