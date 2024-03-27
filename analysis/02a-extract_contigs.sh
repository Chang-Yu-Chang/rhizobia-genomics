#!/usr/bin/env zsh
source ~/.zshrc
source 00-env_vars.sh

# This script splits each genome into contigs

mamba activate samtools
mkdir -p $folder_genomics/contigs

# Assembled genomes
for i in {1..38}; do
    genome_fa=$folder_genomics/genomes/$genome_ids[$i].fasta
    contig_names=($(grep -e "^>" $genome_fa | sed 's/^>//'))

    # Use samtools faidx to extract each contig and save it to a separate file
    for cn in "${contig_names[@]}"; do
        contig_length=$(samtools faidx $genome_fa $cn | wc -c)
        if [[ $contig_length -gt 10000 ]]
            then
                samtools faidx $genome_fa $cn > $folder_genomics/contigs/${genome_ids[$i]}_${cn}.fasta
            fi
    done
done

# Reference genomes
mamba activate seqtk

## Undo soft wrap
cd $folder_genomics/genomes
mkdir -p $folder_genomics/genomes/genomes_ref/
for ref in em1021 wsm419 casidaa; do
    grep '>' $ref.fasta | sed 's/^[^>]*>//g' > genomes_ref/$ref'_contigs'.txt
    seqtk subseq $ref.fasta genomes_ref/$ref'_contigs'.txt > $folder_genomics/genomes/genomes_ref/$ref.fasta
done

# Split the contigs
ref='em1021'
genome_fa=$folder_genomics/genomes/genomes_ref/$ref.fasta
grep -A 1 "NC_003047" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_chromosome'.fasta
grep -A 1 "pSymA" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_psyma'.fasta
grep -A 1 "pSymB" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_psymb'.fasta

ref='wsm419'
genome_fa=$folder_genomics/genomes/genomes_ref/$ref.fasta
grep -A 1 "NC_009636" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_chromosome'.fasta
grep -A 1 "pSMED01" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_psmed01'.fasta
grep -A 1 "pSMED02" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_psmed02'.fasta
grep -A 1 "pSMED03" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_psmed03'.fasta

ref='casidaa'
genome_fa=$folder_genomics/genomes/genomes_ref/$ref.fasta
grep -A 1 "chromosome" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_chromosome'.fasta
grep -A 1 "pCasidaAA" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_pcasidaaa'.fasta
grep -A 1 "pCasidaAB" $genome_fa | tail -n 1 > $folder_genomics/contigs/$ref'_pcasidaab'.fasta
