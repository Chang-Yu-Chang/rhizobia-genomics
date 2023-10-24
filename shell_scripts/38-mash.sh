# Mash: Fast genome and metagenome distance estimation using MinHash

cd
source ~/.zshrc

# Check
mamba activate mash
mamba env list

medaka_consensus=$1
mash_folder=$2
refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"


# Split the consensus into individual contigs
cd $mash_folder
awk '/^>/{if(x) close(x); x=substr($0,2) ".fasta"} {print > x}' $medaka_consensus

# Screening a read set for containment of RefSeq genomes
mash screen -p 10 -w $refseq_db $medaka_consensus  > "$mash_folder/screen.tab"
# `-p 10` Parallelism. This many threads will be spawned for processing
# `-w`  Winner-takes-all strategy for identity estimates. After counting hashes for each query, hashes that appear in multiple queries will be removed from all except the one with the best identity (ties broken by larger query), and other identities will be reduced. This removes output redundancy, providing a rough compositional outline.
## Repeat for contigs
for ct in ./contig*.fasta
do
    mash screen -p 10 -w $refseq_db $ct > "${ct/.fasta/.tab}"
done





# # Sketch the reads, using -m 2 to improve results by ignoring single-copy k-mers, which are more likely to be erroneous:
# mash sketch -m 2 $medaka_consensus -o "$mash_folder/consensus.fasta.msh"
# # `-m 2` Minimum copies of each k-mer required to pass noise filter for reads. Implies -r. [1]
# # `-o <path>`  Output prefix (first input file used if unspecified). The suffix '.msh' will be appended.
# ## Repeat for contigs
# for ct in ./contig*.fasta
# do
#     mash sketch -m 2 $ct -o "$ct.msh"
# done
#
# # Run mash dist with the RefSeq archive as the reference and the read sketch as the query:
# mash dist $refseq_db "$mash_folder/consensus.fasta.msh" > "$mash_folder/distances.tab"
# # `mash dist`
# ## Repeat for contigs
# for ct in ./contig*.msh
# do
#     mash dist $refseq_db $ct > "${ct/.fasta.msh/.tab}"
# done
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
