#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
echo "07-pangenome"





"THIS ONE IS NOT FINISHED YET"


# Calculate ANI
mamba activate fastani

folder_consensus="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus/summary/34-medaka"
#fastANI -q $folder_consensus/consensus_g1.fasta -r $folder_consensus/consensus_g2.fasta -o test.out
realpath $folder_consensus/*.fasta > "$folder_consensus/list_consensus.txt"
fastANI -t 10 \
    --ql "$folder_consensus/list_consensus.txt" \
    --rl "$folder_consensus/list_consensus.txt" \
    -o $folder_consensus/ani.out
# `--ql` list of names of query sequences in fasta
# `--rl` list of names of reference sequences in fasta



# 4. Generate an anvi’o genomes storage
# 4.1 Create a input list of genomes
mkdir "$folder_anvio/04-genomes"
cd "$folder_anvio/04-genomes"
echo -e "name\tcontigs_db_path" > external_genomes.txt
for i in Chang_Q5C_{1..19} usda1106 em1021 em1022 wsm419; do; echo -e "$i\t$folder_anvio/03-contigs_db/$i.db" >> external_genomes.txt; done
cat external_genomes.txt # this should have 19 + 4 = 23 genomes

# 4.2 Generate one storage database for holding the 22 genomes
anvi-gen-genomes-storage \
    -e "$folder_anvio/04-genomes/external_genomes.txt" \
    -o "$folder_anvio/04-genomes/ensifer-GENOMES.db" \
    --gene-caller Prodigal
# -e
# --gene-caller uses the external gene caller included in each genome's db
# -o it must ends with -GENOMES.db



# 5. Run a pangenome analysis
pangenome_id="05-pangenome"
mkdir -p "$folder_anvio/$pangenome_id"

anvi-pan-genome \
    --force-overwrite \
    -g "$folder_anvio/04-genomes/ensifer-GENOMES.db" \
    -n "ensifer" \
    --output-dir "$folder_anvio/$pangenome_id" \
    --num-threads 20 \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --genome-names Chang_Q5C_2,Chang_Q5C_3,Chang_Q5C_4,Chang_Q5C_5,Chang_Q5C_6,Chang_Q5C_8,Chang_Q5C_9,Chang_Q5C_10,Chang_Q5C_11,Chang_Q5C_13,Chang_Q5C_15,Chang_Q5C_16,Chang_Q5C_17,Chang_Q5C_19,usda1106,em1021,em1022,wsm419 \
    --min-occurrence 1
# -g anvio genomes storage file
# -n project name
# --minbit the minimum minibit value
# --mcl-inflation MCL inflation parameter
# --genome-names use a subset of genomes
# --min-occurrence 2 removes singletons
# If you publish results from this workflow, please do not forget to cite DIAMOND (doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL (http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)

# 4.1 add the ANI metric among genomes
anvi-import-misc-data \
    "$folder_data/temp/plasmidsaurus/summary/34-medaka/ani_ensifer.txt" \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    --target-data-table layers \
    --target-data-group ANI_percent_identity
# -e A two-column TAB-delimited flat text file that lists anvi'o contigs databases. Same txt used to make the anvio database
# -p useful when a anivo pangenome is available. The ANI will be written into the pangenome database
# --program program for computing ANI
# Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, please do not forget to properly credit their work.

# 4.3 Check the pangenome database
anvi-db-info "$folder_anvio/$pangenome_id/ensifer-PAN.db"

# 5. Display the pangenome analysis
anvi-display-pan \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    -g "$folder_anvio/04-genomes/ensifer-GENOMES.db"
# -p anvio pan database. The output of anvi-pan-genome
# -g anvio genomes storage file. The output of anvi-gen-genomes-storage

# Then the bin collection has to be saved for anvi-summarize to work


# 6. Explore pangenome
# 6.1 Split the pangenome
# Once you have bins labeled in the pan genome and saved into a "collection" called default
# In this case I have three bins: core, better_core, and duplicated gene pairs
anvi-split \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    -g "$folder_anvio/04-genomes/ensifer-GENOMES.db" \
    -C default \
    -o "$folder_anvio/$pangenome_id/split_pangenome"
# -C collection that contains the bins. the collection of bins has to be specificed in the GUI
# -o output folder

# 6.2 Display the split pangenome
anvi-display-pan \
    -p "$folder_anvio/$pangenome_id/split_pangenome/duplicated_gene_pair/PAN.db" \
    -g "$folder_anvio/04-genomes/ensifer-GENOMES.db"

# 6.3 Provide the trait data as a layer
# The provide trait dataset is only the growth traits and it's from 44_plot_anvio.R
anvi-import-misc-data \
    "$folder_data/temp/42-isolates_anvio.txt" \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    --target-data-table layers

# 6.4 Remove a layer in case it's needed
# anvi-delete-misc-data \
#     -p "$folder_anvio/pangenome/ensifer-PAN.db" \
#     --target-data-table layers \
#     --keys-to-remove r_catogory


# 6.5 compute function enrichment
# This is not yet possible becasue the T
# anvi-compute-functional-enrichment-in-pan \
#     -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
#     -g "$folder_anvio/ensifer-GENOMES.db" \
#     --category r_category \
#     --annotation-source Prokka:Prodigal \
#     -o "$folder_anvio/enriched_functions.txt"
# This only works for categorical traits
# This uses a Rscript to calculate the enrichment score. Need to make sure the R packages are included
# It seems that my mamba environment has it, but the Rscript CLT uses the base Rscript instead of the mamba env bin/Rscript
# Do this to install the required R packages in base for now
# install.packages("tidyverse")
# install.packages("optparse")
# BiocManager::install("qvalue")


# 7. Summarize the pangenomes
anvi-summarize \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    -g "$folder_anvio/04-genomes/ensifer-GENOMES.db" \
    -o "$folder_anvio/$pangenome_id/summary"
#    -C default \
# -C the collection of bins has to be specificed in the GUI using anvi-display-pan

# unzip the txt file
gunzip "$folder_anvio/$pangenome_id/summary/ensifer_gene_clusters_summary.txt.gz"
