#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script QCs the assemblies

cd $folder_shell

for i in {1..32}
do
    echo $genome_ids[$i]
    dir=$folder_genomics/assembly/$genome_ids[$i]

    # Quality control the assemblies
    mamba activate quast
    quast \
        $dir/medaka/consensus.fasta \
        -t 10 \
        -o $folder_genomics/assembly/$genome_ids[$i]/quast

done





