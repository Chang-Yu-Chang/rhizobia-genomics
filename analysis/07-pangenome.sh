#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script implements pangenome analysis

cd $folder_shell
echo "07-pangenome"
mkdir -p $folder_genomics/pangenome
folder_pangenome="$folder_genomics/pangenome"

# Roary
rm -rf $folder_pangenome/roary_gff
mkdir -p $folder_pangenome/roary_gff/
for i in {1..23} {24..41}
do
    cp "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka/annotated.gff" "$folder_pangenome/roary_gff/$sample_ids[$i].gff"
done

zsh 07f-roary.sh \
    "$folder_pangenome/roary_rhizobia" \
    "$folder_pangenome/roary_gff"

# Panaroo
mkdir -p $folder_pangenome/panaroo_gff/

for i in {2..6} {8..11} 13 {15..17}, {19..37} {38..42}
do
    cp "$folder_genomes/$sample_ids[$i]/05-gene_annotation/prokka/annotated.gff" "$folder_pangenome/panaroo_gff/$sample_ids[$i].gff"
done


mamba activate panaroo
panaroo -i $folder_pangenome/panaroo_gff/*.gff \
    -o $folder_pangenome/panaroo_ensifer \
    -t 10 \
    --clean-mode moderate --remove-invalid-genes






# 4.3 Check the pangenome database
#anvi-db-info "$folder_anvio/$pangenome_id/ensifer-PAN.db"

# 5. Display the pangenome analysis
# anvi-display-pan \
#     -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
#     -g "$folder_anvio/04-genomes/ensifer-GENOMES.db"
# -p anvio pan database. The output of anvi-pan-genome
# -g anvio genomes storage file. The output of anvi-gen-genomes-storage

# Then the bin collection has to be saved for anvi-summarize to work


# 6. Explore pangenome
# 6.1 Split the pangenome
# Once you have bins labeled in the pan genome and saved into a "collection" called default
# In this case I have three bins: core, better_core, and duplicated gene pairs
# anvi-split \
#     -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
#     -g "$folder_anvio/04-genomes/ensifer-GENOMES.db" \
#     -C default \
#     -o "$folder_anvio/$pangenome_id/split_pangenome"
# -C collection that contains the bins. the collection of bins has to be specificed in the GUI
# -o output folder

# 6.2 Display the split pangenome
# anvi-display-pan \
#     -p "$folder_anvio/$pangenome_id/split_pangenome/duplicated_gene_pair/PAN.db" \
#     -g "$folder_anvio/04-genomes/ensifer-GENOMES.db"

# 6.3 Provide the trait data as a layer
# The provide trait dataset is only the growth traits and it's from 44_plot_anvio.R
# anvi-import-misc-data \
#     "$folder_data/temp/42-isolates_anvio.txt" \
#     -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
#     --target-data-table layers

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




# # Generate an anvi’o genomes storage
# # Create a input list of external genomes
# echo -e "name\tcontigs_db_path" > "$folder_pangenome/list_genomes_for_anvio.txt"
# for i in {1..41}
# do
#     echo "$sample_ids[$i]\t$folder_genomes/$sample_ids[$i]/06-pangenome_prep/genome.db" >> "$folder_pangenome/list_genomes_for_anvio.txt"
# done
#
# # Generate one storage database for holding the 23 genomes
# zsh 07b-create_anvio_pan.sh \
#     "$folder_pangenome/list_genomes_for_anvio.txt" \
#     "$folder_pangenome/rhizobia-GENOMES.db"
#
#
# # Run a pangenome analysis
# # For all rhizobia genomes
# for i in {1..41}
# do
#     echo $sample_ids[$i]
# done >| "$folder_pangenome/list_genomes_rhizobia.txt"
#
# zsh 07c-run_anvio_pan.sh \
#     "$folder_pangenome/rhizobia-GENOMES.db" \
#     "rhizobia" \
#     "$folder_pangenome/rhizobia" \
#     "$folder_pangenome/list_genomes_rhizobia.txt"
#
# # For only ensifer spp
# for i in {2..6} {8..11} 13 {15..17} {19..23} {24..41}
# do
#     echo $sample_ids[$i]
# done >| "$folder_pangenome/list_genomes_ensifer.txt"
#
# zsh 07c-run_anvio_pan.sh \
#     "$folder_pangenome/rhizobia-GENOMES.db" \
#     "ensifer" \
#     "$folder_pangenome/ensifer" \
#     "$folder_pangenome/list_genomes_ensifer.txt"
#
#
# # Calculate ANI
# for i in {1..23} {24..41}
# do
#     echo -e "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"
# done >| "$folder_pangenome/list_genomes_for_ani_rhizobia.txt"
#
# zsh 07a-fastani.sh \
#     "$folder_pangenome/list_genomes_for_ani_rhizobia.txt" \
#     "$folder_pangenome/fastani_rhizobia.txt"
#
# for i in {2..6} {8..11} 13 {15..17} {19..23} {24..41}
# do
#     echo -e "$folder_genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta"
# done >| "$folder_pangenome/list_genomes_for_ani_ensifer.txt"
#
# zsh 07a-fastani.sh \
#     "$folder_pangenome/list_genomes_for_ani_ensifer.txt" \
#     "$folder_pangenome/fastani_ensifer.txt"
#
# # Add the ANI metric among genomes
# zsh 07d-import_ani.sh \
#     "$folder_pangenome/fastani_rhizobia.out" \
#     "$folder_pangenome/rhizobia/rhizobia-PAN.db"
#
# zsh 07d-import_ani.sh \
#     "$folder_pangenome/fastani_ensifer.out" \
#     "$folder_pangenome/ensifer/ensifer-PAN.db"
#
# # Summarize the pangenomes
# zsh 07e-summarize_pan.sh \
#     "$folder_pangenome/rhizobia/rhizobia-PAN.db" \
#     "$folder_pangenome/rhizobia-GENOMES.db" \
#     "$folder_pangenome/rhizobia/summary" \
#     "rhizobia"
#
# zsh 07e-summarize_pan.sh \
#     "$folder_pangenome/ensifer/ensifer-PAN.db" \
#     "$folder_pangenome/rhizobia-GENOMES.db" \
#     "$folder_pangenome/ensifer/summary" \
#     "ensifer"