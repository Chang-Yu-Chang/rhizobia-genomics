#!/usr/bin/env zshs
source ~/.zshrc

# This estimates the genome distance to refseq
# $1: genome in fasta
# $2: mash_folder
# $3: refseq_db

conda activate
mamba activate mash

# Split the consensus into individual contigs
cd $2
awk '/^>/{if(x) close(x); x=substr($0,2) ".fasta"} {print > x}' $1

# Screening a read set for containment of RefSeq genomes
mash screen -p 10 -w $3 $1  > "$2/screen.tab"
# `-p 10` Parallelism. This many threads will be spawned for processing
# `-w`  Winner-takes-all strategy for identity estimates. After counting hashes for each query, hashes that appear in multiple queries will be removed from all except the one with the best identity (ties broken by larger query), and other identities will be reduced. This removes output redundancy, providing a rough compositional outline.

## Repeat for contigs
for ct in $2/c_*.fasta
do
    mash screen -p 10 -w $3 $ct > "${ct/.fasta/.tab}"
done
