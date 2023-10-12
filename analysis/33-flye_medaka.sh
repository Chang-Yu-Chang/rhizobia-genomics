cd
source ~/.zshrc


# Using information acquired from the Miniasm assembly, re-downsample the reads to ~100x coverage (do nothing if there isn't at least 100x coverage) with heavy weight applied to removing low quality reads (helps small plasmids stick around)
mamba activate filtlong
mamba env list

filtered_reads=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filtlong/filter_reads.fastq.gz
downsample_reads2=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka/downsample_reads.fatq
assemble_fq=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/miniasm_pilot.fasta
genome_size=$(grep -v ">" "$assemble_fq" | grep -E -o "G|C|T|A|N" | wc -l)

filtlong --target_bases $((genome_size * 100)) --mean_q_weight 10 -a "$assemble_fq" "$filtered_reads" | > "$downsample_reads2"




# Select for high quality ONT reads
mamba activate flye
mamba env list

#assemble_paf=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/miniasm_pilot.paf.gz
#assemble_gfa=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/miniasm_pilot.gfa
assemble_fq=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/02-miniasm/miniasm_pilot.fasta
assembly_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/03-flye_medaka

flye --nano-raw "$assemble_fq" --out-dir "$assembly_folder" --threads 4
