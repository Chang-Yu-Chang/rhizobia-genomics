cd
source ~/.zshrc


# Annotate genomes via bakta
mamba activate bakta
mamba env list

bakta_database=$1
medaka_consensus=$2
bakta_folder=$3

# bakta_database=~/bioinformatics/bakta/db # This database is mandatory and must be downloaded before annotation
# medaka_consensus=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/04-medaka/consensus.fasta
# bakta_folder=~/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/05-bakta

# bakta --db $bakta_database --verbose --output $bakta_folder --thread 8 $medaka_consensus
# `--db "$bakta_database"` database directory
# `--verbose` Print verbose information
# `--output "$bakta_folder"` output directory
# `--thread 8`
# `"$medaka_consensus"`  Genome sequences in (zipped) fasta format


# Show the circos plot
bakta_plot --output $bakta_folder --prefix circos "$bakta_folder/consensus.json"

#bakta_plot --output test --prefix test --config config.yaml --sequences 1,2 input.json
