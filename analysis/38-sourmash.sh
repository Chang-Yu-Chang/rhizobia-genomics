cd
source ~/.zshrc


# Check
mamba activate sourmash
mamba env list

medaka_consensus=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.fasta
sourmash_dir=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/08-sourmash/
sourmash_sig=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/08-sourmash/consensus.fasta.sig
genbank_db=~/bioinformatics/sourmash/genbank-k31.lca.json.gz

# Create a signature for this genome
sourmash sketch dna --param-string scaled=1000,k=31 "$medaka_consensus" --output-dir "$sourmash_dir"
# `sketch dna` command reads in DNA sequences and outpouts DNA sketches
# `--param-string scaled=1000,k=31` signature parameters to use
# `--output-dir` output computed signatures to this directory

# Classify the signature with sourmash lca classify
sourmash lca classify --db "$genbank_db" --query "$sourmash_sig"
# `lca classify` classify genomes
# `--db` databases to use to classify
# `--query` query signatures to classify

# You can also summarize the taxonomic distribution of the content with lca summarize:
sourmash lca summarize --db "$genbank_db" --query "$sourmash_sig"
# `lca summarize` summarize mixture
