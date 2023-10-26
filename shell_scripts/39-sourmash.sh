cd
source ~/.zshrc


# Check
mamba activate sourmash
mamba env list

medaka_consensus=$1
sourmash_folder=$2
# medaka_consensus="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_1/04-medaka/consensus.fasta"
# sourmash_folder="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_1/09-sourmash"
gtdb_db="/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip"
sourmash_sig="$sourmash_folder/consensus.sig"
gathered_csv="$sourmash_folder/gathered.csv"

# Create a signature for this genome
sourmash sketch dna --check-sequence -f -p scaled=1000,k=31 $medaka_consensus -o "$sourmash_folder/consensus.sig"

# `sketch dna` command reads in DNA sequences and outpouts DNA sketches
# `-p scaled=1000,k=31` signature parameters to use
# `-o` output computed signatures to this directory

# We recommend using the Zipfile databases for sourmash gather and the SBT databases for sourmash search.
sourmash gather $sourmash_sig $gtdb_db -o $gathered_csv

# Classify the signature
#sourmash tax --gather-csv HSMA33MX_gather_x_gtdbrs202_k31.csv --taxonomy gtdb-rs202.taxonomy.v2.csv
#sourmash tax genome --gather-csv $gather_csv --taxonomy "$sourmash_folder/taxonomy.csv"



# Select the best reference genomes to use
# Evaluate containment, that is, what fraction of the read content is contained in the genome:
#sourmash search -k 31 $sourmash_sig $genbank_db
# Selet the reference genomes to use for a metagenome analysis
#sourmash gather $sourmash_sig $genbank_db
# `gather` selects the best reference genomes to use for a metagenome analysis, by finding the smallest set of non-overlapping matches to the query in a database.  This is specifically meant for metagenome and genome bin analysis.
# it takes two positional arguments `query` and `databases`

















# sourmash lca classify --db $genbank_db --query $sourmash_sig
# # `lca classify` classify genomes
# # `--db` databases to use to classify
# # `--query` query signatures to classify

# You can also summarize the taxonomic distribution of the content with lca summarize:
#sourmash lca summarize --db $genbank_db --query $sourmash_sig
# `lca summarize` summarize mixture




