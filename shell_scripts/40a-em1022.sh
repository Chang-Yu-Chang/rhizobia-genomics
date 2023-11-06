#
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
mamba activate ncbi-datasets

# Search NCBI for, EM1021 and EM1022, WSM419
#' EM1021 was SM1021 or 1021; EM1022 was WSM1022
datasets summary genome taxon 'ensifer meliloti' --assembly-source refseq --as-json-lines | \
dataformat tsv genome --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name | \
grep '1021'
# GCF_000006965.1

datasets summary genome taxon 'ensifer meliloti' --assembly-source refseq --as-json-lines | \
dataformat tsv genome --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name | \
grep '1022'
# GCF_013315775.1

datasets summary genome taxon 'ensifer medicae' --assembly-source refseq --as-json-lines | \
dataformat tsv genome --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name | \
grep 'WSM419'
# GCF_000017145.1


# Download the EM1021 and EM1022 genomes
datasets download genome accession GCF_000006965.1 GCF_013315775.1 GCF_000017145.1 \
--include genome,cds,gff3,seq-report \
--filename  "$folder_data/raw/ensifer_ncbi.zip"
cd "$folder_data/raw"
rm README.md
unzip ensifer_ncbi.zip
mv ncbi_dataset ensifer_ncbi

# Move the files
mkdir -p "$folder_data/temp/ncbi"
cd "$folder_data/temp/ncbi"
mkdir -p "em1021/prokka"
mkdir -p "em1022/prokka"
mkdir -p "wsm419/prokka"

cp "$folder_data/raw/ensifer_ncbi/data/GCF_000006965.1/GCF_000006965.1_ASM696v1_genomic.fna" \
"$folder_data/temp/ncbi/em1021/genome.fasta"
cp "$folder_data/raw/ensifer_ncbi/data/GCF_013315775.1/GCF_013315775.1_ASM1331577v1_genomic.fna" \
"$folder_data/temp/ncbi/em1022/genome.fasta"
cp "$folder_data/raw/ensifer_ncbi/data/GCF_000017145.1/GCF_000017145.1_ASM1714v1_genomic.fna" \
"$folder_data/temp/ncbi/wsm419/genome.fasta"


# Reformate ncbi genomes to contigs compatible to anvio
mamba activate anvio-8
anvi-script-reformat-fasta "$folder_data/temp/ncbi/em1021/genome.fasta" -o "$folder_data/temp/ncbi/em1021/contigs.fasta" --simplify-names
anvi-script-reformat-fasta "$folder_data/temp/ncbi/em1022/genome.fasta" -o "$folder_data/temp/ncbi/em1022/contigs.fasta" --simplify-names
anvi-script-reformat-fasta "$folder_data/temp/ncbi/wsm419/genome.fasta" -o "$folder_data/temp/ncbi/wsm419/contigs.fasta" --simplify-names


# Annotate these genomes using prokka
mamba activate prokka
prokka --force --kingdom Bacteria --prefix annotated --gcode 11 \
--outdir "$folder_data/temp/ncbi/em1021/prokka" "$folder_data/temp/ncbi/em1021/contigs.fasta"

prokka --force --kingdom Bacteria --prefix annotated --gcode 11 \
--outdir "$folder_data/temp/ncbi/em1022/prokka" "$folder_data/temp/ncbi/em1022/contigs.fasta"

prokka --force --kingdom Bacteria --prefix annotated --gcode 11 \
--outdir "$folder_data/temp/ncbi/wsm419/prokka" "$folder_data/temp/ncbi/wsm419/contigs.fasta"



# Pan genome analysis thes twoe using roary -> in the next script
















