cd
source ~/.zshrc

conda activate

folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus"
mkdir -p "$folder_temp_result/summary/04-medaka"


# Copy
for i in {1..19}
do
    cp "$folder_temp_result/Chang_Q5C_$i/04-medaka/consensus.fasta" \
    "$folder_temp_result/summary/04-medaka/consensus_g$i.fasta"
done
