#!/usr/bin/env zsh
source ~/.zshrc

# This prepares the anvio database for pangenome analysis
# Following this tutorial https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome
# $1: genome in fasta
# $2: mash_folder
# $3: refseq_db

conda activate
mamba activate anvio-8
# Check setup
#anvi-self-test --suite pangenomics # You might have to install any dependency according to the error message if any

# folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
# mkdir -p "$folder_data/temp/anvio"
# folder_anvio="$folder_data/temp/anvio"
#cd $folder_anvio

for i in {1..23}
do


done

# 2.2 Importing prokka annotations into anvio. Prokka should have had happened
mamba activate anvio-8
cd ~/bioinformatics/anvio # For using the gff_parser.py
mkdir -p "$folder_anvio/02-prokka/gene_calls"
mkdir -p "$folder_anvio/02-prokka/gene_annot"

for i in Chang_Q5C_{1..19} usda1106 em1021 em1022 wsm419
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

for i in Chang_Q5C_{1..19} usda1106 em1021 em1022 wsm419
do
    # Create the contig database attached with external gene call
    anvi-gen-contigs-database \
        --force-overwrite \
        -f "$folder_anvio/01-fasta/$i.fasta" \
        -o "$folder_anvio/03-contigs_db/$i.db" \
        --external-gene-calls "$folder_anvio/02-prokka/gene_calls/$i.txt" \
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
#anvi-db-info "$folder_anvio/03-contigs_db/Chang_Q5C_15.db"

# 2.3 Run hmms, which helps annotate the genes in your contigs-db. This adds the HMM info to the database
for i in Chang_Q5C_{1..19} usda1106 em1021 em1022 wsm419; do; anvi-run-hmms -c "$folder_anvio/03-contigs_db/$i.db"; done

# 2.4 output the contig stat from each contigs-db
mkdir -p "$folder_anvio/03-contigs_db/contigs_stat"
for i in Chang_Q5C_{1..19} usda1106 em1021 em1022 wsm419
do
    anvi-display-contigs-stats \
        "$folder_anvio/03-contigs_db/$i.db" \
        --report-as-text \
        -o "$folder_anvio/03-contigs_db/contigs_stat/$i.txt"
        # --report-as-text  If you give this flag, Anvi'o will not open new browser to show Contigs database statistics and write all stats to TAB separated file and you should also give --output-file with this flag otherwise Anvi'o will complain.
done

# # 2.5 Export contigs
# mkdir -p "$folder_anvio/contigs"
# for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
# do
#     anvi-export-contigs \
#         -c "$folder_anvio/genomes/$i.db" \
#         -o "$folder_anvio/contigs/$i.fasta"
# done



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
