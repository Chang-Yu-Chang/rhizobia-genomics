cd
source ~/.zshrc

cd ~/Desktop/lab/local-adaptation/shell_scripts

conda activate

#batch_id="Chang_Q5C_results_repeated"
# batch_id="Chang_Q5C_results"
# sample_id="Chang_Q5C_1"
# zsh 00c-script_for_one.sh $batch_id $sample_id

# Nanocompare all raw reads
mamba activate nanocomp
zsh 30a-nanocomp.sh

# Run the master script
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_1"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_2"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_3"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_4"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_5"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_6"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_7"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_8"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_9"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_10"
zsh 00c-script_for_one.sh "Chang_Q5C_results_repeated" "Chang_Q5C_11"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_12"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_13"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_14"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_15"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_16"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_17"
zsh 00c-script_for_one.sh "Chang_Q5C_results_repeated" "Chang_Q5C_18"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_19"

# Extract the read length and ASC text
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_1"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_2"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_3"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_4"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_5"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_6"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_7"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_8"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_9"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_10"
zsh 31a-extract_reads.sh "Chang_Q5C_results_repeated" "Chang_Q5C_11"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_12"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_13"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_14"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_15"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_16"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_17"
zsh 31a-extract_reads.sh "Chang_Q5C_results_repeated" "Chang_Q5C_18"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_19"


# After all analysis above
# Copy assembled genomes to one folder
zsh 34a-move_consensus.sh
# Make busco plot
zsh 37a-plot_busco.sh
# Pangenome analysis
zsh 41-roary.sh
zsh 42-scoary.sh



















