#!/usr/bin/env zsh
source ~/.zshrc

# This script searched and downloads the NCBI genomes

conda activate
mamba activate ncbi-datasets

# Read the mapping file
table_file=$1
batch_names=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 1)}")
sample_ids=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 2)}")
accessions=("${(@f)$(tail -n +2 $table_file | cut -d ',' -f 3)}")

# Download the Ensifer genomes
datasets download genome accession $accessions[38,41] \
    --include genome,gff3,gbff \
    --filename "$2/ensifer_ncbi.zip"
# GCF_002197065.1 is E. meliloti usda1106
# GCF_000006965.1 is E. meliloti em1021
# GCF_013315775.1 is E. meliloti em1022
# GCF_000017145.1 is E. medicae wsm419

# Clean up the files and names
cd $2
unzip ensifer_ncbi.zip
rm -rf ncbi_genomes
mv ncbi_dataset/data ncbi_genomes
rm -rf README.md ensifer_ncbi.zip ncbi_dataset


# Move the genome files
mamba activate anvio-8

for i in {38..41}
do
    # Move files
    mkdir -p "$2/genomes/$sample_ids[$i]/02-denovo_assembly/ncbi"
    cp $2/ncbi_genomes/$accessions[$i]/*_genomic.fna $2/genomes/$sample_ids[$i]/02-denovo_assembly/ncbi/genome.fasta
    cp $2/ncbi_genomes/$accessions[$i]/genomic.gbff $2/genomes/$sample_ids[$i]/02-denovo_assembly/ncbi/genome.gbff

    # Reformate ncbi genomes to contigs compatible to anvio
    anvi-script-reformat-fasta \
        "$2/genomes/$sample_ids[$i]/02-denovo_assembly/ncbi/genome.fasta" \
        -o "$2/genomes/$sample_ids[$i]/02-denovo_assembly/genome.fasta" \
        --simplify-names \
        --report-file "$2/genomes/$sample_ids[$i]/02-denovo_assembly/cleaned_names.txt" \
        --min-len 500
    # --simplify-names edit deflines to make sure they have simple names. This is important for BAM
    # --report-file reports the changes to deflines
    # --min-len 500 minimum length of contigs to keep
done



# Also move ncbi strains data to a medaka folder for calculating ani
#for i in em1021 em1022 wsm419; do; cp "$folder_data/temp/plasmidsaurus/$i/contigs.fasta" "$folder_data/temp/plasmidsaurus/summary/34-medaka/$i.fasta"; done
# # Search NCBI for, EM1021, EM1022, USDA1106, WSM419
# # EM1021 was SM1021 or 1021; EM1022 was WSM1022
# datasets summary genome taxon 'ensifer meliloti' --assembly-source refseq --as-json-lines | \
# dataformat tsv genome --fields accession,annotinfo-release-date,organism-name,organism-infraspecific-strain | \
# grep -e "USDA1106" -e "WSM1022" -e "1021"
# # GCF_002197065.1 is usda1105
# # GCF_000006965.1 is em1021
# # GCF_013315775.1 is em1022
#
# datasets summary genome taxon 'ensifer medicae' --assembly-source refseq --as-json-lines | \
# dataformat tsv genome --fields accession,annotinfo-release-date,organism-name | \
# grep -e 'WSM419'
# # GCF_000017145.1 is wsm419
