cd
source ~/.zshrc


# Check
mamba activate mash
mamba env list

medaka_consensus=$1
mash_folder=$2

refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"
mash_sketch="$mash_folder/consensus.mash"
mash_distance="$mash_folder/distances.tab"


# medaka_consensus=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.fasta
# mash_sketch=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/07-mash/consensus.msh
# refseq_db=~/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh
# mash_distance=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/07-mash/distances.tab

# Split the consensus into individual contigs
cd $mash_folder
awk '/^>/{if(x) close(x); x=substr($0,2) ".fasta"} {print > x}' $medaka_consensus

# Sketch the reads, using -m 2 to improve results by ignoring single-copy k-mers, which are more likely to be erroneous:
for ct in ./contig*.fasta
do
    mash sketch -m 2 $ct -o "${ct/.fasta/.msh}"
done
# `-m 2` Minimum copies of each k-mer required to pass noise filter for reads. Implies -r. [1]
# `-o <path>`  Output prefix (first input file used if unspecified). The suffix '.msh' will be appended.

# Run mash dist with the RefSeq archive as the reference and the read sketch as the query:
for ct in ./contig*.msh
do
    mash dist -s $refseq_db -p 0.1 $ct > "${ct/.msh/.tab}"
done
#mash dist -s $refseq_db -p 0.1 "$mash_sketch.msh" > $mash_distance
# `mash dist`
# `>`

# Sort the results to see the top hits and their p-values:
#sort -gk3 $mash_distance | head
# `-g` sort by general numerical values
# `-k3` sort omitting the first and second fields

# NOT DONE YET. I NEED TO FIGURE HOW TO FIND THE TAXONOMY


