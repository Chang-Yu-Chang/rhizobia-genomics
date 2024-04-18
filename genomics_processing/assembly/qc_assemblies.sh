#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script QCs the assemblies

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

    # Checkm
    mamba activate checkm
    checkm lineage_wf \
        -x fasta \
        -t 10 \
        $dir/medaka/ \
        $folder_genomics/assembly/$genome_ids[$i]/checkm
done

# Checkm
# mamba activate checkm
# checkm lineage_wf \
#     -t 10 \
#     -x fasta \
#     $folder_genomics/fasta/genomes/ \
#     $folder_genomics/fasta/genomes/checkm
