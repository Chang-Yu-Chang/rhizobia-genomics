cd
source ~/.zshrc


# Check
mamba activate mash
mamba env list

medaka_consensus=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.fasta
refseq_db=~/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh
mash_sketch=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/07-mash/consensus.fasta
mash_distance=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/07-mash/distances.tab

# Sketch the reads, using -m 2 to improve results by ignoring single-copy k-mers, which are more likely to be erroneous:
mash sketch -m 2 "$medaka_consensus" -o "$mash_sketch"
# `-m 2` Minimum copies of each k-mer required to pass noise filter for reads. Implies -r. [1]
# `-o <path>`  Output prefix (first input file used if unspecified). The suffix '.msh' will be appended.

# Run mash dist with the RefSeq archive as the reference and the read sketch as the query:
mash dist "$refseq_db" "$mash_sketch.msh" > "$mash_distance"
# `mash dist`
# `>`

# Sort the results to see the top hits and their p-values:
sort -gk3 "$mash_distance" | head
# `g`


# NOT DONE YET. I NEED TO FIGURE HOW TO FIND THE TAXONOMY
