cd
source ~/.zshrc

conda activate
mamba activate nanocomp
mamba env list

folder_raw_result="/Users/cychang/Dropbox/lab/local-adaptation/data/raw"

NanoComp --verbose \
--fastq "$folder_raw_result/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_2/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_3/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_4/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_5/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_6/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_7/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_8/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_9/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_10/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results_repeated/Chang_Q5C_11/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_12/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_13/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_14/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_15/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_16/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_17/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results_repeated/Chang_Q5C_18/reads/raw_reads.fastq.gz" \
"$folder_raw_result/Chang_Q5C_results/Chang_Q5C_19/reads/raw_reads.fastq.gz" \
--names g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15 g16 g17 g18 g19 \
--outdir "/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus"
