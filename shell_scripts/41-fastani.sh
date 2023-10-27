cd
source ~/.zshrc

# Check
mamba activate fastani
mamba env list


#
folder_consensus="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/04-medaka"

fastANI -q $folder_consensus/consensus_g1.fasta -r $folder_consensus/consensus_g2.fasta -o test.out

realpath $folder_consensus/* > list_consensus.txt

fastANI --ql list_consensus.txt --rl list_consensus.txt -o $folder_consensus/test.out
