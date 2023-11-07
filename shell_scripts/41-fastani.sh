cd
source ~/.zshrc

mamba activate fastani
mamba env list

folder_consensus="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/34-medaka"
#fastANI -q $folder_consensus/consensus_g1.fasta -r $folder_consensus/consensus_g2.fasta -o test.out
realpath $folder_consensus/*.fasta > "$folder_consensus/list_consensus.txt"
fastANI -t 10 \
    --ql "$folder_consensus/list_consensus.txt" \
    --rl "$folder_consensus/list_consensus.txt" \
    -o $folder_consensus/ani.out
# `--ql` list of names of query sequences in fasta
# `--rl` list of names of reference sequences in fasta
