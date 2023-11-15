#!/usr/bin/env zsh
# Following this tutorial https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome


folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
mkdir -p "$folder_data/temp/anvio"
folder_anvio="$folder_data/temp/anvio"
cd $folder_anvio

# 0. Check setup
mamba activate anvio-8
#anvi-self-test --suite pangenomics # You might have to install any dependency according to the error message if any

# 1. Clean the assembled genome fasta to have simplified deflines
mkdir -p "$folder_anvio/01-fasta"

for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    # Clean the fasta
    anvi-script-reformat-fasta \
        "$folder_data/temp/plasmidsaurus/$i/04-medaka/consensus.fasta" \
        -o "$folder_anvio/01-fasta/$i.fasta" \
        --simplify-names \
        --report-file "$folder_anvio/01-fasta/$i.txt" \
        --min-len 500
done

# --simplify-names edit deflines to make sure they have simple names. This is important for BAM
# --report-file reports the changes to deflines
# --min-len 500 minimum length of contigs to keep


# 2. Prokka
# 2.1 Run Prokka
mamba activate prokka
mkdir -p "$folder_anvio/02-prokka"

for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    prokka \
        --force --kingdom Bacteria --prefix annotated --gcode 11 \
        --outdir "$folder_anvio/02-prokka/$i" \
        "$folder_anvio/01-fasta/$i.fasta"
done
# --force Force overwriting existing output folder (default OFF)
# --kingdom [X]      Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
#  --prefix [X]       Filename output prefix [auto] (default '')
# --gcode [N]        Genetic code / Translation table (set if --kingdom is set) (default '0')

# 2.2 Importing prokka annotations into anvio. Prokka should have had happened
mamba activate anvio-8
cd ~/bioinformatics/anvio # For using the gff_parser.py
mkdir -p "$folder_anvio/02-prokka/gene_calls"
mkdir -p "$folder_anvio/02-prokka/gene_annot"

for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    python3 gff_parser.py \
        "$folder_anvio/02-prokka/$i/annotated.gff" \
        --gene-calls "$folder_anvio/02-prokka/gene_calls/$i.txt" \
        --annotation "$folder_anvio/02-prokka/gene_annot/$i.txt"
done

# 3. Create one contigs database per assembly FASTA file
mamba activate anvio-8
# 3.1 Create the contigs-db
mkdir -p "$folder_anvio/03-contigs_db"

for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    # Create the contig database attached with external gene call
    anvi-gen-contigs-database \
        --force-overwrite \
        -f "$folder_anvio/01-fasta/$i.fasta" \
        -o "$folder_anvio/03-contigs_db/$i.db" \
        --external-gene-calls "$folder_anvio/gene_calls/$i.txt" \
        --project-name $i \
        --ignore-internal-stop-codons
    # --external-gene-calls A TAB-delimited file to define external gene calls.
    # --ignore-internal-stop-codons This is only relevant when you have an external gene calls file

    # Import the functional annotations
    anvi-import-functions \
        -c "$folder_anvio/03-contigs_db/$i.db" \
        -i "$folder_anvio/02-prokka/gene_annot/$i.txt"
done


# 2.2 Examine a anvio database if needed
anvi-db-info "$folder_anvio/03-contigs_db/Chang_Q5C_15.db"

# 2.3 Run hmms, which helps annotate the genes in your contigs-db. This adds the HMM info to the database
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419; do; anvi-run-hmms -c "$folder_anvio/03-contigs_db/$i.db"; done

# 2.4 output the contig stat from each contigs-db
mkdir -p "$folder_anvio/contigs_stat"
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    anvi-display-contigs-stats \
        "$folder_anvio/genomes/$i.db" \
        --report-as-text \
        -o "$folder_anvio/contigs_stat/$i.txt"
done

# # 2.5 Export contigs
# mkdir -p "$folder_anvio/contigs"
# for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
# do
#     anvi-export-contigs \
#         -c "$folder_anvio/genomes/$i.db" \
#         -o "$folder_anvio/contigs/$i.fasta"
# done









# 3. Generate an anvi’o genomes storage
# 3.1 Create a input list of genomes
cd "$folder_anvio"
echo -e "name\tcontigs_db_path" > external_genomes.txt
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419; do; echo -e "$i\t$folder_anvio/genomes/$i.db" >> external_genomes.txt; done
cat external_genomes.txt # this should have 19 + 3 = 22 genomes

# 3.2 Generate one storage database for holding the 22 genomes
anvi-gen-genomes-storage \
    -e "$folder_anvio/external_genomes.txt" \
    -o "$folder_anvio/ensifer-GENOMES.db" \
    --gene-caller Prodigal
# `-e`
# `--gene-caller` uses the external gene caller included in each genome's db
# `-o` it must ends with -GENOMES.db



# 4. Run a pangenome analysis
pangenome_id="pangenome"
mkdir -p "$folder_anvio/$pangenome_id"
anvi-pan-genome \
    --force-overwrite \
    -g "ensifer-GENOMES.db" \
    -n "ensifer" \
    --output-dir "$folder_anvio/$pangenome_id" \
    --num-threads 20 \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --genome-names Chang_Q5C_2,Chang_Q5C_3,Chang_Q5C_4,Chang_Q5C_5,Chang_Q5C_6,Chang_Q5C_8,Chang_Q5C_9,Chang_Q5C_10,Chang_Q5C_11,Chang_Q5C_13,Chang_Q5C_15,Chang_Q5C_16,Chang_Q5C_17,Chang_Q5C_19,em1021,em1022,wsm419 \
    --min-occurrence 1
# `-g` anvio genomes storage file
# `-n` project name
# `--minbit` the minimum minibit value
# `--mcl-inflation` MCL inflation parameter
# `--genome-names` use a subset of genomes
# `--min-occurrence 2` removes singletons
# If you publish results from this workflow, please do not forget to cite DIAMOND (doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL (http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)

# 4.1 add the ANI metric among genomes
anvi-import-misc-data "$folder_data/temp/plasmidsaurus/summary/34-medaka/ani_ensifer.txt" \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    --target-data-table layers \
    --target-data-group ANI_percent_identity
# `-e` A two-column TAB-delimited flat text file that lists anvi'o contigs databases. Same txt used to make the anvio database
# `-p` useful when a anivo pangenome is available. The ANI will be written into the pangenome database
# `--program` program for computing ANI
# Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, please do not forget to properly credit their work.

# 4.3 Check the pangenome database
anvi-db-info "$folder_anvio/$pangenome_id/ensifer-PAN.db"


# 5. Display the pangenome analysis
anvi-display-pan \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    -g "$folder_anvio/ensifer-GENOMES.db"
# `-p` anvio pan database. The output of anvi-pan-genome
# `-g` anvio genomes storage file. The output of anvi-gen-genomes-storage

# Then the bin collection has to be saved for anvi-summarize to work


# 6. Explore pangenome
# 6.1 Split the pangenome
# Once you have bins labeled in the pan genome and saved into a "collection" called default
# In this case I have three bins: core, better_core, and duplicated gene pairs
anvi-split \
    -p "$folder_anvio/$pangenome_id/ensifer-PAN.db" \
    -g "$folder_anvio/ensifer-GENOMES.db" \
    -C default \
    -o "$folder_anvio/split_$pangenome_id"
# `-C` collection that contains the bins. the collection of bins has to be specificed in the GUI
# `-o` output folder

# 6.2 Desplay the split pangenome
anvi-display-pan \
    -p "$folder_anvio/split_$pangenome_id/duplicated_gene_pair/PAN.db" \
    -g "$folder_anvio/ensifer-GENOMES.db"

# 6.3 Provide the trait data as a layer
# The provide trait dataset is only the growth traits and it's from 44_plot_anvio.R
anvi-import-misc-data "$folder_data/temp/42-isolates_anvio.txt" \
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
    -g "$folder_anvio/ensifer-GENOMES.db" \
    -C default \
    -o "$folder_anvio/summary_$pangenome_id"
# `-C` the collection of bins has to be specificed in the GUI using anvi-display-pan

# unzip the txt file
gunzip "$folder_anvio/summary_$pangenome_id/ensifer_gene_clusters_summary.txt.gz"

# convert svg to pdf
# folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
# cd "$folder_anvio/figures"
# rsvg-convert -f pdf -o ensifer.pdf ensifer.svg
# rsvg-convert -f pdf -o ensifer_circle.pdf ensifer_circle.svg
# rsvg-convert -f pdf -o ensifer_duplicated.pdf ensifer_duplicated.svg

