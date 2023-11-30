cd
source ~/.zshrc

cd ~/Desktop/lab/local-adaptation/shell_scripts
source 00-env_vars.sh
conda activate

#batch_id="Chang_Q5C_results_repeated"
# batch_id="Chang_Q5C_results"
# sample_id="Chang_Q5C_1"
# zsh 00c-script_for_one.sh $batch_id $sample_id

# Nanocompare all raw reads
# mamba activate nanocomp
# zsh 30a-nanocomp.sh

# Download NCBI Ensifer genomes
datasets download genome accession GCF_002197065.1 GCF_000006965.1 GCF_013315775.1 GCF_000017145.1 \
    --include genome,gff3,gbff \
    --filename  "$folder_data/raw/ensifer_ncbi.zip"
cd "$folder_data/raw"
unzip ensifer_ncbi.zip
rm README.md

# Move the files
for i in usda1106 em1021 em1022 wsm419
do
    mkdir -p "$folder_genomes/$i/05-gene_annotation"
    mkdir -p "$folder_genomes/$i/06-pangenome_prep"
done

cp $folder_data/raw/ncbi_dataset/data/GCF_002197065.1/*_genomic.fna $folder_genomes/usda1106/genome.fasta
cp $folder_data/raw/ncbi_dataset/data/GCF_000006965.1/*_genomic.fna $folder_genomes/em1021/genome.fasta
cp $folder_data/raw/ncbi_dataset/data/GCF_013315775.1/*_genomic.fna $folder_genomes/em1022/genome.fasta
cp $folder_data/raw/ncbi_dataset/data/GCF_000017145.1/*_genomic.fna $folder_genomes/wsm419/genome.fasta

cp $folder_data/raw/ncbi_dataset/data/GCF_002197065.1/genomic.gbff $folder_genomes/usda1106/genome.gbff
cp $folder_data/raw/ncbi_dataset/data/GCF_000006965.1/genomic.gbff $folder_genomes/em1021/genome.gbff
cp $folder_data/raw/ncbi_dataset/data/GCF_013315775.1/genomic.gbff $folder_genomes/em1022/genome.gbff
cp $folder_data/raw/ncbi_dataset/data/GCF_000017145.1/genomic.gbff $folder_genomes/wsm419/genome.gbff

# Clean the genome fasta contig names to have simplified deflines
mamba activate anvio-8
for i in usda1106 em1021 em1022 wsm419
do
    anvi-script-reformat-fasta \
        "$folder_genomes/$i/genome.fasta" \
        -o "$folder_genomes/$i/06-pangenome_prep/genome.fasta" \
        --simplify-names \
        --report-file "$folder_genomes/$i/06-pangenome_prep/report.txt" \
        --min-len 500
done

# Index the Reference Genome. For long reads, tools like Minimap2 or NGMLR are often preferred
mamba activate minimap2
mkdir -p "$folder_genomics/reference"
for i in usda1106 em1021 em1022 wsm419
do
    minimap2 -d "$folder_genomics/reference/$i.mmi" "$folder_genomes/$i/06-pangenome_prep/genome.fasta"
done



# Run the master script
for i in Chang_Q5C_{1..10} Chang_Q5C_{11..17} Chang_Q5C_19
do
    zsh 00c-script_for_one.sh "Chang_Q5C_results" $i
done

for i in Chang_Q5C_11 Chang_Q5C_18
do
    zsh 00c-script_for_one.sh "Chang_Q5C_results_repeated" $i
done





zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_1"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_2"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_3"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_4"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_5"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_6"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_7"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_8"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_9"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_10"
zsh 00c-script_for_one.sh "Chang_Q5C_results_repeated" "Chang_Q5C_11"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_12"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_13"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_14"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_15"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_16"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_17"
zsh 00c-script_for_one.sh "Chang_Q5C_results_repeated" "Chang_Q5C_18"
zsh 00c-script_for_one.sh "Chang_Q5C_results" "Chang_Q5C_19"

# Extract the read length and ASC text
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_1"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_2"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_3"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_4"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_5"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_6"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_7"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_8"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_9"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_10"
zsh 31a-extract_reads.sh "Chang_Q5C_results_repeated" "Chang_Q5C_11"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_12"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_13"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_14"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_15"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_16"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_17"
zsh 31a-extract_reads.sh "Chang_Q5C_results_repeated" "Chang_Q5C_18"
zsh 31a-extract_reads.sh "Chang_Q5C_results" "Chang_Q5C_19"


# After all analysis above
# Copy assembled genomes to one folder
zsh 34a-move_consensus.sh
# Make busco plot
zsh 37a-plot_busco.sh
# Pangenome analysis
zsh 41-roary.sh
zsh 42-scoary.sh



















