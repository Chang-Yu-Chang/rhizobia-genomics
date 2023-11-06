# Following this tutorial https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome


folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
mkdir -p "$folder_data/temp/anvio"
cd "$folder_data/temp/anvio"

# Check setup
mamba activate anvio-8
anvi-self-test --suite pangenomics # You might have to install any dependency according to the error message if any

# Reformate ncbi genomes to contigs compatible to anvio
anvi-script-reformat-fasta "$folder_data/temp/ncbi/em1021/genome.fasta" -o "$folder_data/temp/ncbi/em1021/contigs.fasta" --simplify-names
anvi-script-reformat-fasta "$folder_data/temp/ncbi/em1022/genome.fasta" -o "$folder_data/temp/ncbi/em1022/contigs.fasta" --simplify-names
anvi-script-reformat-fasta "$folder_data/temp/ncbi/wsm419/genome.fasta" -o "$folder_data/temp/ncbi/wsm419/contigs.fasta" --simplify-names

# Importing prokka annotations into anvio. Prokka should have taken place
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

# Create a contigs database from a FASTA file
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
        --external-gene-calls "$folder_data/temp/anvio/$i/anvio/gene_calls.txt" \
        --name "$i" \
        --ignore-internal-stop-codons
    # Import the functional annotations
    anvi-import-functions \
        -c "$folder_data/temp/anvio/genomes/$i.db" \
        -i "$folder_data/temp/ncbi/$i/anvio/gene_annot.txt"
done


# Create a input list of genomes
cd "$folder_data/temp/anvio"
echo -e "name\tcontigs_db_path" > external_genomes.txt
for i in {1..19}; do; echo -e "Chang_Q5C_$i\t$folder_data/temp/anvio/genomes/Chang_Q5C_$i.db" >> external_genomes.txt; done
for i in em1021 em1022 wsm419; do; echo -e "$i\t$folder_data/temp/anvio/genomes/$i.db" >> external_genomes.txt; done
cat external_genomes.txt

# Generate an anvi’o genomes storage
anvi-gen-genomes-storage -e external_genomes.txt -o ensifer_genomes.db

# Run a pangenome analysis
anvi-pan-genome -g ensifer_genomes.db -n "ensifer_pangenome" \
    --output-dir "$folder_data/temp/anvio" \
    --num-threads 6 \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --use-ncbi-blast
# `-g` anvio genomes storage file
# `-n` project name
# `--minbit` the minimum minibit value
# `--mcl-inflation` MCL inflation parameter
# `--genome-names` use a subset of genomes


# Display the pangenome analysis
anvi-display-pan -p PROJECT-PAN.db \
    -g PROJECT-PAN-GENOMES.db \
    --title

# ``














