#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script performs reference-guided genome assembly

cd $folder_shell
echo "02-guided_assembly"

# Download NCBI genomes
zsh 02a-download_ncbi_genomes.sh

#

for i in {1..19}
do
    zsh 02a-miniasm.sh  \
        "$folder_genomics/reference/usda1106.mmi" \
        $sample_id


    zsh 02b-align_genomes.sh  \
        "$folder_genomics/reference/usda1106.mmi" \
        $sample_id

done


