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
    # mamba activate quast
    # quast \
    #     $dir/medaka/consensus.fasta \
    #     -t 10 \
    #     -o $folder_genomics/assembly/$genome_ids[$i]/quast

    # BUSCO
    # mamba activate busco
    # busco \
    #     -i $dir/medaka/consensus.fasta \
    #     -m genome \
    #     -c 10 \
    #     -l alphaproteobacteria_odb10 \
    #     -o "" \
    #     --out_path $folder_genomics/assembly/$genome_ids[$i]/busco

    # Define BUSCO lineages
    mamba activate busco
    lineages=(
        #"alphaproteobacteria_odb10"
        "rhizobiales_odb10"
        "rhizobiaceae_odb12"
        "sinorhizobium_odb12"
    )

    for lineage in $lineages
    do
        echo "  -> BUSCO lineage: $lineage"

        if [ -d $dir/busco/$lineage ]; then
            echo "  -> Skipping $lineage (results already exist)"
            continue
        fi

        busco \
            -i $dir/medaka/consensus.fasta \
            -m genome \
            -c 10 \
            -l $lineage \
            --out_path $dir/busco/$lineage \
            -o ""
    done
done
