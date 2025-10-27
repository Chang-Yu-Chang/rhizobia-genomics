#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script assesses the assembly quality using quast, busco, and checkm

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

# Calculate percentage removed
mamba activate seqkit
for i in {1..38}
do
    id=${genome_ids[$i]}
    dir=$folder_genomics/assembly/$id
    infile=$dir/medaka/consensus.fasta
    outfile=$folder_genomics/fasta/$id.fasta

    # Total genome length before filtering
    total_before=$(seqkit stats -T $infile | awk 'NR==2 {print $5}')
    # Filter contigs shorter than 100 kb
    seqkit seq -m 100000 -g $infile > $outfile
    # Total genome length after filtering
    total_after=$(seqkit stats -T $outfile | awk 'NR==2 {print $5}')
    # Calculate % removed
    percent_removed=$(awk -v before=$total_before -v after=$total_after 'BEGIN {printf("%.2f", (before - after)/before * 100)}')
    echo "$id: $(printf "%'.0f" $total_before) bp → $(printf "%'.0f" $total_after) bp  (Removed: $percent_removed%)"
done


# Checkm
mamba activate checkm
checkm lineage_wf \
    -t 8 \
    -x fasta \
    $folder_genomics/fasta/ \
    $folder_genomics/checkm

# Extract summary table
checkm qa \
  $folder_genomics/checkm/lineage.ms \
  $folder_genomics/checkm/ \
  -o 2 -t 8 > $folder_genomics/checkm/checkm_summary.tsv
