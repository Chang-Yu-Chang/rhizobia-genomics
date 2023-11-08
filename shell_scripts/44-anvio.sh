# Following this tutorial https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome


folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
mkdir -p "$folder_data/temp/anvio"
cd "$folder_data/temp/anvio"

# 0. Check setup
mamba activate anvio-8
anvi-self-test --suite pangenomics # You might have to install any dependency according to the error message if any


# 1. Importing prokka annotations into anvio. Prokka should have had happened
cd ~/bioinformatics/anvio
for i in {1..19}
do
    mkdir -p "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/11-anvio"
    python3 gff_parser.py "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/10-prokka/annotated.gff" \
    --gene-calls "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/11-anvio/gene_calls.txt" \
    --annotation "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/11-anvio/gene_annot.txt"
done

for i in em1021 em1022 wsm419
do
    mkdir -p "$folder_data/temp/ncbi/$i/anvio"
    python3 gff_parser.py "$folder_data/temp/ncbi/$i/prokka/annotated.gff" \
        --gene-calls "$folder_data/temp/ncbi/$i/anvio/gene_calls.txt" \
        --annotation "$folder_data/temp/ncbi/$i/anvio/gene_annot.txt"
done

# 2. Create one contigs database per assembly FASTA file
mkdir -p "$folder_data/temp/anvio/genomes"
for i in {1..19}
do
    # Create the contig database attached with external gene call
    anvi-gen-contigs-database --force-overwrite \
        -f "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/04-medaka/consensus.fasta" \
        -o "$folder_data/temp/anvio/genomes/Chang_Q5C_$i.db" \
        --external-gene-calls "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/11-anvio/gene_calls.txt" \
        --project-name "Chang_Q5C_$i" \
        --ignore-internal-stop-codons
    # Import the functional annotations
    anvi-import-functions \
        -c "$folder_data/temp/anvio/genomes/Chang_Q5C_$i.db" \
        -i "$folder_data/temp/plasmidsaurus/Chang_Q5C_$i/11-anvio/gene_annot.txt"

done

for i in em1021 em1022 wsm419
do
    # Create the contig database attached with external gene call
    anvi-gen-contigs-database --force-overwrite \
        -f "$folder_data/temp/ncbi/$i/contigs.fasta" \
        -o "$folder_data/temp/anvio/genomes/$i.db" \
        --external-gene-calls "$folder_data/temp/ncbi/$i/anvio/gene_calls.txt" \
        --project-name "$i" \
        --ignore-internal-stop-codons
    # Import the functional annotations
    anvi-import-functions \
        -c "$folder_data/temp/anvio/genomes/$i.db" \
        -i "$folder_data/temp/ncbi/$i/anvio/gene_annot.txt"
done

# 2.1 Examine a anvio database if needed
# i=1
# anvi-db-info "$folder_data/temp/anvio/genomes/Chang_Q5C_$i.db"

# 2.2 Run hmms, which helps annotate the genes in your contigs-db. This add the HMM info to the database
for i in {1..19}; do; anvi-run-hmms -c "$folder_data/temp/anvio/genomes/Chang_Q5C_$i.db"; done
for i in em1021 em1022 wsm419; do; anvi-run-hmms -c "$folder_data/temp/anvio/genomes/$i.db"; done


# 3. Generate an anvi’o genomes storage
# 3.1 Create a input list of genomes
cd "$folder_data/temp/anvio"
echo -e "name\tcontigs_db_path" > external_genomes.txt
for i in {1..19}; do; echo -e "Chang_Q5C_$i\t$folder_data/temp/anvio/genomes/Chang_Q5C_$i.db" >> external_genomes.txt; done
for i in em1021 em1022 wsm419; do; echo -e "$i\t$folder_data/temp/anvio/genomes/$i.db" >> external_genomes.txt; done
cat external_genomes.txt # this should have 19 + 3 = 22 genomes

# 3.2 Generate one storage file for holding the 22 genomes
anvi-gen-genomes-storage \
    -e "$folder_data/temp/anvio/external_genomes.txt" \
    -o "$folder_data/temp/anvio/ensifer-GENOMES.db" \
    --gene-caller Prodigal
# `-e`
# `--gene-caller` uses the external gene caller included in each genome's db
# `-o` it must ends with -GENOMES.db


# 4. Run a pangenome analysis
mkdir -p "$folder_data/temp/anvio/pangenome"
anvi-pan-genome \
    --force-overwrite \
    -g ensifer-GENOMES.db \
    -n "ensifer" \
    --output-dir "$folder_data/temp/anvio/pangenome" \
    --num-threads 20 \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --genome-names Chang_Q5C_2,Chang_Q5C_3,Chang_Q5C_4,Chang_Q5C_5,Chang_Q5C_6,Chang_Q5C_8,Chang_Q5C_9,Chang_Q5C_10,Chang_Q5C_11,Chang_Q5C_13,Chang_Q5C_15,Chang_Q5C_16,Chang_Q5C_17,Chang_Q5C_19,em1021,em1022,wsm419 \
    --min-occurrence 2
#    --use-ncbi-blast
# `-g` anvio genomes storage file
# `-n` project name
# `--minbit` the minimum minibit value
# `--mcl-inflation` MCL inflation parameter
# `--genome-names` use a subset of genomes
# `--min-occurrence 2` removes singletons
# If you publish results from this workflow, please do not forget to cite DIAMOND (doi:10.1038/nmeth.3176), unless you use it with --use-ncbi-blast flag, and MCL (http://micans.org/mcl/ and doi:10.1007/978-1-61779-361-5_15)

# 4.1 add the ANI metric among genomes
anvi-import-misc-data "$folder_data/temp/plasmidsaurus/summary/34-medaka/ani_ensifer.txt" \
    -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
    --target-data-table layers \
    --target-data-group ANI_percent_identity
# anvi-compute-genome-similarity \
#     -e "$folder_data/temp/anvio/external_genomes.txt" \
#     -o "$folder_data/temp/anvio/ani" \
#     -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
#     --program pyANI
# # `-e` A two-column TAB-delimited flat text file that lists anvi'o contigs databases. Same txt used to make the anvio database
# # `-p` useful when a anivo pangenome is available. The ANI will be written into the pangenome database
# # `--program` program for computing ANI
# # Anvi'o will use 'PyANI' by Pritchard et al. (DOI: 10.1039/C5AY02550H) to compute ANI. If you publish your findings, please do not forget to properly credit their work.

# Check the pangenome database
anvi-db-info "$folder_data/temp/anvio/pangenome/ensifer-PAN.db"


# 5. Display the pangenome analysis
anvi-display-pan \
    -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
    -g "$folder_data/temp/anvio/ensifer-GENOMES.db"
# `-p` anvio pan database. The output of anvi-pan-genome
# `-g` anvio genomes storage file. The output of anvi-gen-genomes-storage


# 6. Explore pangenome
# 6.1 Split the pangenome
# Once you have bins labeled in the pan genome and saved into a "collection" called default
# In this case I have three bins: core, better_core, and duplicated gene pairs
anvi-split \
    -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
    -g "$folder_data/temp/anvio/ensifer-GENOMES.db" \
    -C default \
    -o "$folder_data/temp/anvio/split_pangenome"
# `-C` collection that contains the bins
# `-o` output folder
anvi-display-pan \
    -p "$folder_data/temp/anvio/split_pangenome/duplicated_gene_pair/PAN.db" \
    -g "$folder_data/temp/anvio/ensifer-GENOMES.db"



#anvi-display-pan -p "$folder_data/temp/anvio/split_pangenome/core/PAN.db" -g "$folder_data/temp/anvio/ensifer-GENOMES.db"

# 6.2 Provide the trait data as a layer
anvi-import-misc-data "$folder_data/temp/42-isolates_anvio.txt" \
    -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
    --target-data-table layers
# The provide trait dataset is only the growth traits and it's from 44_plot_anvio.R
# Rmove layer
# anvi-delete-misc-data \
#     -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
#     --target-data-table layers \
#     --keys-to-remove r_catogory


# 6.3 compute function enrichment
anvi-compute-functional-enrichment-in-pan \
    -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
    -g "$folder_data/temp/anvio/ensifer-GENOMES.db" \
    --category r_category \
    --annotation-source Prokka:Prodigal \
    -o "$folder_data/temp/anvio/enriched_functions.txt"
# This only works for categorical traits
# This uses a Rscript to calculate the enrichment score. Need to make sure the R packages are included
# It seems that my mamba environment has it, but the Rscript CLT uses the base Rscript instead of the mamba env bin/Rscript
# Do this to install the required R packages in base for now
# install.packages("tidyverse")
# install.packages("optparse")
# BiocManager::install("qvalue")


# 7. Summarize the pangenomes
anvi-summarize \
    -p "$folder_data/temp/anvio/pangenome/ensifer-PAN.db" \
    -g "$folder_data/temp/anvio/ensifer-GENOMES.db" \
    -C default \
    -o "$folder_data/temp/anvio/summary"


# unzip the txt file
gunzip "$folder_data/temp/anvio/summary/ensifer_gene_clusters_summary.txt.gz"

# convert svg to pdf
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
cd "$folder_data/temp/anvio/figures"
rsvg-convert -f pdf -o ensifer.pdf ensifer.svg
rsvg-convert -f pdf -o ensifer_circle.pdf ensifer_circle.svg
rsvg-convert -f pdf -o ensifer_duplicated.pdf ensifer_duplicated.svg









