#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script splits each genome into contigs

mamba activate samtools
mkdir -p $folder_genomics/fasta/contigs

## Assembled genomes
for i in {1..32}; do
    genome_fa=$folder_genomics/genomes/$genome_ids[$i].fasta
    contig_names=($(grep -e "^>" $genome_fa | sed 's/^>//'))

    # Use samtools faidx to extract each contig and save it to a separate file
    for cn in "${contig_names[@]}"; do
        contig_length=$(samtools faidx $genome_fa $cn | wc -c)
        # Remove contigs shorter than 10k
        if [[ $contig_length -gt 10000 ]]
            then
                samtools faidx $genome_fa $cn > $folder_genomics/fasta/contigs/${genome_ids[$i]}_${cn}.fasta
            fi
    done
done

## Reference genomes
mamba activate seqtk

## Undo soft wrap
mkdir -p $folder_genomics/fasta/genomes/genomes_ref/
cd $folder_genomics/fasta/genomes
for ref in em1021 em1022 usda1106 wsm419 casidaa; do
    grep '>' $ref.fasta | sed 's/^[^>]*>//g' > genomes_ref/$ref'_contigs'.txt
    seqtk subseq $ref.fasta genomes_ref/$ref'_contigs'.txt > genomes_ref/$ref.fasta
done

# Split the contigs
ref='em1021'
genome_fa=$folder_genomics/fasta/genomes/genomes_ref/$ref.fasta
grep -A 1 "NC_003047" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_chromosome'.fasta
grep -A 1 "pSymA" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psyma'.fasta
grep -A 1 "pSymB" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psymb'.fasta

ref='em1022'
genome_fa=$folder_genomics/fasta/genomes/genomes_ref/$ref.fasta
grep -A 1 "chromosome" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_chromosome'.fasta
grep -A 1 "pA" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psyma'.fasta
grep -A 1 "pB" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psymb'.fasta

ref='usda1106'
genome_fa=$folder_genomics/fasta/genomes/genomes_ref/$ref.fasta
grep -A 1 "chromosome" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_chromosome'.fasta
grep -A 1 "psymA" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psyma'.fasta
grep -A 1 "psymB" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psymb'.fasta

ref='wsm419'
genome_fa=$folder_genomics/fasta/genomes/genomes_ref/$ref.fasta
grep -A 1 "NC_009636" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_chromosome'.fasta
grep -A 1 "pSMED01" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psmed01'.fasta
grep -A 1 "pSMED02" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psmed02'.fasta
grep -A 1 "pSMED03" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_psmed03'.fasta

ref='casidaa'
genome_fa=$folder_genomics/fasta/genomes/genomes_ref/$ref.fasta
grep -A 1 "chromosome" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_chromosome'.fasta
grep -A 1 "pCasidaAA" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_pcasidaaa'.fasta
grep -A 1 "pCasidaAB" $genome_fa | tail -n 2 > $folder_genomics/fasta/contigs/$ref'_pcasidaab'.fasta

