cd
source ~/.zshrc


# Check
mamba activate sourmash
mamba env list

medaka_consensus=$1
sourmash_folder=$2
genbank_db=$3
sourmash_sig=$4
gather_csv=$5

# medaka_consensus=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.fasta
# sourmash_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/08-sourmash
# sourmash_sig=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/08-sourmash/consensus.fasta.sig
# genbank_db=~/bioinformatics/sourmash/genbank-k31.lca.json.gz

# Create a signature for this genome
sourmash sketch dna --param-string scaled=1000,k=31 $medaka_consensus --output-dir $sourmash_folder
# `sketch dna` command reads in DNA sequences and outpouts DNA sketches
# `--param-string scaled=1000,k=31` signature parameters to use
# `--output-dir` output computed signatures to this directory

# Select the best reference genomes to use
sourmash gather $sourmash_sig $genbank_db -o $gather_csv "THIS IS NOT FINISHED YET"
# it takes two positional arguments `query` and `databases`
# `-o` output CSV containg matches to this file

# Classify the signature with sourmash
sourmash tax --gather-csv HSMA33MX_gather_x_gtdbrs202_k31.csv --taxonomy gtdb-rs202.taxonomy.v2.csv
sourmash tax genome -sourma$medaka_consensus --gather-csv test_gat.csv --taxonomy test_tax.csv

# sourmash lca classify --db $genbank_db --query $sourmash_sig
# # `lca classify` classify genomes
# # `--db` databases to use to classify
# # `--query` query signatures to classify

# You can also summarize the taxonomic distribution of the content with lca summarize:
sourmash lca summarize --db $genbank_db --query $sourmash_sig
# `lca summarize` summarize mixture


