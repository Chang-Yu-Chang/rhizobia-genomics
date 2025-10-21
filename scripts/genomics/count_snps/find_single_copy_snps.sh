#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# Single copy core genes
mamba activate biopython

mkdir -p $folder_data/genomics_analysis/count_snps
for i in 4 5 6 8 9 10 11 13 16 17 19 {20..27} 29 30 31 32 33 34 35 36 37 39 40 41 42 43 44 45;
do
    echo g$i
done | > $folder_data/genomics_analysis/count_snps/list_genomes_symbiotic.txt

# Count snps
python find_single_copy_snps.py \
    $folder_genomics/pangenome/aligned_gene_sequences/ \
    $folder_data/genomics_analysis/count_snps/list_genomes_symbiotic.txt \
    $folder_data/genomics_analysis/count_snps/list_snps_symbiotic.txt

# mamba activate seqkit
# seqkit stats rpmI.aln.fas
# seqkit stats rplW.aln.fas
