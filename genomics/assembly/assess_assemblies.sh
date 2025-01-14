#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assesses the quality of assembly

cd $folder_shell

for i in {1..38}
do
    echo $genome_ids[$i]
    dir=$folder_genomics/assembly/$genome_ids[$i]

    # Quast
    mamba activate quast
    quast \
        $dir/medaka/consensus.fasta \
        -t 10 \
        -o $folder_genomics/assembly/$genome_ids[$i]/quast

    # BUSCO
    mamba activate busco
    busco \
        -i $dir/medaka/consensus.fasta \
        -m genome \
        -l alphaproteobacteria_odb10 \
        -o "" \
        --out_path $folder_genomics/assembly/$genome_ids[$i]/busco

done
