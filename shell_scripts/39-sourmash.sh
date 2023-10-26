cd
source ~/.zshrc


# Check
mamba activate sourmash
mamba env list

medaka_consensus=$1
sourmash_folder=$2
# medaka_consensus="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_1/04-medaka/consensus.fasta"
# sourmash_folder="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_1/09-sourmash"
genbank_db="/Users/cychang/bioinformatics/sourmash/genbank-2022.03-bacteria-k31.zip"
#"THIS SHOULD BE A SBT FILE"
sourmash_sig="$sourmash_folder/consensus.fasta.sig"
gather_csv="$sourmash_folder/gather.csv"

# medaka_consensus="~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_1/04-medaka/consensus.fasta"
# sourmash_folder="~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_1/09-sourmash"

# Create a signature for this genome
sourmash sketch dna --param-string scaled=1000,k=31 $medaka_consensus --output-dir $sourmash_folder
# `sketch dna` command reads in DNA sequences and outpouts DNA sketches
# `--param-string scaled=1000,k=31` signature parameters to use
# `--output-dir` output computed signatures to this directory

# Selet the reference genomes to use for a metagenome analysis
#sourmash gather $sourmash_sig $genbank_db
# `gather` selects the best reference genomes to use for a metagenome analysis, by finding the smallest set of non-overlapping matches to the query in a database.  This is specifically meant for metagenome and genome bin analysis.
# it takes two positional arguments `query` and `databases`

# We recommend using the Zipfile databases for sourmash gather and the SBT databases for sourmash search.
sourmash search $sourmash_sig $genbank_db --containment

# Classify the signature with sourmash
#sourmash tax --gather-csv HSMA33MX_gather_x_gtdbrs202_k31.csv --taxonomy gtdb-rs202.taxonomy.v2.csv
sourmash tax genome --gather-csv $gather_csv --taxonomy "$sourmash_folder/taxonomy.csv"



# Select the best reference genomes to use
# Evaluate containment, that is, what fraction of the read content is contained in the genome:
#sourmash search -k 31 $sourmash_sig $genbank_db

















# sourmash lca classify --db $genbank_db --query $sourmash_sig
# # `lca classify` classify genomes
# # `--db` databases to use to classify
# # `--query` query signatures to classify

# You can also summarize the taxonomic distribution of the content with lca summarize:
#sourmash lca summarize --db $genbank_db --query $sourmash_sig
# `lca summarize` summarize mixture




