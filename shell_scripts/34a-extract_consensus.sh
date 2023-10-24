cd
source ~/.zshrc

mamba activate bioawk
#mamba env list

# Raw reads
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
folder_raw_result="/Users/cychang/Dropbox/lab/local-adaptation/data/raw/$1"
folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/$1"
sample_id=$2

sample_id=$2

echo $sample_id

# consensus_fas="$folder_temp_result/$sample_id/04-medaka/consensus.fasta"
# consensus_tab="$folder_temp_result/$sample_id/04-medaka/consensus.txt"
consensus_fas="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.fasta"
consensus_tab="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.txt"

bioawk -c fastx '{print $name, $seq}' $consensus_fas > $consensus_tab
