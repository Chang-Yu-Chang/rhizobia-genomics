folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
mamba activate ncbi-datasets

# 1. Search NCBI for, EM1021, EM1022, USDA1106, WSM419
#' EM1021 was SM1021 or 1021; EM1022 was WSM1022
datasets summary genome taxon 'ensifer meliloti' --assembly-source refseq --as-json-lines | \
dataformat tsv genome --fields accession,annotinfo-release-date,organism-name,organism-infraspecific-strain | \
grep -e "USDA1106" -e "WSM1022" -e "1021"
# GCF_002197065.1 is usda1105
# GCF_000006965.1 is em1021
# GCF_013315775.1 is em1022

datasets summary genome taxon 'ensifer medicae' --assembly-source refseq --as-json-lines | \
dataformat tsv genome --fields accession,annotinfo-release-date,organism-name | \
grep -e 'WSM419'
# GCF_000017145.1 is wsm419



# 2. Download the Ensifer genomes
datasets download genome accession GCF_002197065.1 GCF_000006965.1 GCF_013315775.1 GCF_000017145.1 \
    --include genome,gff3,gbff \
    --filename  "$folder_data/raw/ensifer_ncbi.zip"
cd "$folder_data/raw"
unzip ensifer_ncbi.zip
rm README.md

# Move the files
mkdir -p "$folder_data/temp/plasmidsaurus"
cd "$folder_data/temp/plasmidsaurus"

for i in usda1106 em1021 em1022 wsm419
do
    mkdir -p "$folder_data/temp/plasmidsaurus/$i/10-prokka"
    mkdir -p "$folder_data/temp/plasmidsaurus/$i/04-medaka"
done

cp $folder_data/raw/ncbi_dataset/data/GCF_002197065.1/*_genomic.fna $folder_data/temp/plasmidsaurus/usda1106/genome.fasta
cp $folder_data/raw/ncbi_dataset/data/GCF_000006965.1/*_genomic.fna $folder_data/temp/plasmidsaurus/em1021/genome.fasta
cp $folder_data/raw/ncbi_dataset/data/GCF_013315775.1/*_genomic.fna $folder_data/temp/plasmidsaurus/em1022/genome.fasta
cp $folder_data/raw/ncbi_dataset/data/GCF_000017145.1/*_genomic.fna $folder_data/temp/plasmidsaurus/wsm419/genome.fasta

cp $folder_data/raw/ncbi_dataset/data/GCF_002197065.1/genomic.gbff $folder_data/temp/plasmidsaurus/usda1106/genome.gbff
cp $folder_data/raw/ncbi_dataset/data/GCF_000006965.1/genomic.gbff $folder_data/temp/plasmidsaurus/em1021/genome.gbff
cp $folder_data/raw/ncbi_dataset/data/GCF_013315775.1/genomic.gbff $folder_data/temp/plasmidsaurus/em1022/genome.gbff
cp $folder_data/raw/ncbi_dataset/data/GCF_000017145.1/genomic.gbff $folder_data/temp/plasmidsaurus/wsm419/genome.gbff

cp $folder_data/temp/plasmidsaurus/usda1106/genome.fasta $folder_data/temp/plasmidsaurus/usda1106/04-medaka/consensus.fasta
cp $folder_data/temp/plasmidsaurus/em1021/genome.fasta $folder_data/temp/plasmidsaurus/em1021/04-medaka/consensus.fasta
cp $folder_data/temp/plasmidsaurus/em1022/genome.fasta $folder_data/temp/plasmidsaurus/em1022/04-medaka/consensus.fasta
cp $folder_data/temp/plasmidsaurus/wsm419/genome.fasta $folder_data/temp/plasmidsaurus/wsm419/04-medaka/consensus.fasta




# Also move ncbi strains data to a medaka folder for calculating ani
#for i in em1021 em1022 wsm419; do; cp "$folder_data/temp/plasmidsaurus/$i/contigs.fasta" "$folder_data/temp/plasmidsaurus/summary/34-medaka/$i.fasta"; done

# # Reformate ncbi genomes to contigs compatible to anvio
# mamba activate anvio-8
# anvi-script-reformat-fasta "$folder_data/temp/plasmidsaurus/em1021/genome.fasta" -o "$folder_data/temp/plasmidsaurus/em1021/04-medaka/consensus.fasta" --simplify-names
# anvi-script-reformat-fasta "$folder_data/temp/plasmidsaurus/em1022/genome.fasta" -o "$folder_data/temp/plasmidsaurus/em1022/04-medaka/consensus.fasta" --simplify-names
# anvi-script-reformat-fasta "$folder_data/temp/plasmidsaurus/wsm419/genome.fasta" -o "$folder_data/temp/plasmidsaurus/wsm419/04-medaka/consensus.fasta" --simplify-names
#
# # Annotate these genomes using prokka
# mamba activate prokka
# prokka --force --kingdom Bacteria --prefix annotated --gcode 11 \
# --outdir "$folder_data/temp/plasmidsaurus/em1021/10-prokka" "$folder_data/temp/plasmidsaurus/em1021/04-medaka/consensus.fasta"
#
# prokka --force --kingdom Bacteria --prefix annotated --gcode 11 \
# --outdir "$folder_data/temp/plasmidsaurus/em1022/10-prokka" "$folder_data/temp/plasmidsaurus/em1022/04-medaka/consensus.fasta"
#
# prokka --force --kingdom Bacteria --prefix annotated --gcode 11 \
# --outdir "$folder_data/temp/plasmidsaurus/wsm419/10-prokka" "$folder_data/temp/plasmidsaurus/wsm419/04-medaka/consensus.fasta"
#














# Pan genome analysis thes twoe using roary -> in the next script
