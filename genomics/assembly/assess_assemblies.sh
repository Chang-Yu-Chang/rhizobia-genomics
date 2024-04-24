#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assesses the quality of assembly

cd $folder_shell

for i in {1..32}
do
    echo $genome_ids[$i]
    dir=$folder_genomics/assembly/$genome_ids[$i]

    # Quast
    # mamba activate quast
    # quast \
    #     $dir/medaka/consensus.fasta \
    #     -t 10 \
    #     -o $folder_genomics/assembly/$genome_ids[$i]/quast

    # BUSCO
    mamba activate busco
    busco \
        -i $dir/medaka/consensus.fasta \
        -m genome \
        -l alphaproteobacteria_odb10 \
        -o "" \
        --out_path $folder_genomics/assembly/$genome_ids[$i]/busco

    # Checkm
    # mamba activate checkm
    # checkm lineage_wf \
    #     -x fasta \
    #     -t 10 \
    #     $dir/medaka/ \
    #     $folder_genomics/assembly/$genome_ids[$i]/checkm
done

# Checkm
# mamba activate checkm
# checkm lineage_wf \
#     -t 10 \
#     -x fasta \
#     $folder_genomics/fasta/genomes/ \
#     $folder_genomics/fasta/genomes/checkm


# # BUSCO
# mamba activate busco
# busco \
#     -i $folder_genomics/fasta/genomes/*fasta \
#     -m genome \
#     -l alphaproteobacteria_odb10 \
#     --download_path $folder_genomics/fasta/genomes/busco \
#     -o $genome_ids[$i] \
#     --out_path $folder_genomics/fasta/genomes/busco















