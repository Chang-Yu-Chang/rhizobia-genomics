#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assesses the assembly quality using quast, busco, and checkm

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
    # lineages=(
    #     #"alphaproteobacteria_odb10"
    #     "rhizobiales_odb10"
    #     "rhizobiaceae_odb12"
    #     "sinorhizobium_odb12"
    # )
    #
    # for lineage in $lineages
    # do
    #     echo "  -> BUSCO lineage: $lineage"
    #
    #     if [ -d $dir/busco/$lineage ]; then
    #         echo "  -> Skipping $lineage (results already exist)"
    #         continue
    #     fi
    #
    #     busco \
    #         -i $dir/medaka/consensus.fasta \
    #         -m genome \
    #         -c 10 \
    #         -l $lineage \
    #         --out_path $dir/busco/$lineage \
    #         -o ""
    # done

    # Checkm
    mamba activate checkm
    mkdir -p $dir/checkm/

    # skip if already done
    if [ -d $dir/checkm/ ]; then
        echo "   ⚠️  Skipping (checkm folder exists)"
        continue
    fi

    checkm lineage_wf \
        -t 8 \
        -x fasta \
        $dir/medaka/consensus.fasta \
        $dir/checkm/

done
