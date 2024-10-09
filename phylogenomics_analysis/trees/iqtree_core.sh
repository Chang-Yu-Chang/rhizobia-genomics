#!/usr/bin/env zsh
source ~/.zshrc
source ../../genomics/env_vars.sh

# Create the isolate gene tree
mkdir -p $folder_data/phylogenomics_analysis/trees/mltree/seq_core

# Run IQ-TREE on the combined alignment file
mamba activate iqtree

# iqtree \
#     -nt AUTO \
#     -B 1000 \
#     -s $folder_genomics/pangenome/isolates/core_gene_alignment.aln \
#     --prefix $folder_data/phylogenomics_analysis/trees/mltree/seq_core/seq_core

# -nt AUTO will tell IQ-TREE to automatically determine the best number of cores given the current data and computer.
# -B Replicates for ultrafast bootstrap (>=1000)
elev_med/aligned_gene_sequences/.aln.fas

set_name="elev_med"
for file in $folder_genomics/pangenome/$set_name/aligned_gene_sequences/*.aln.fas; do

    # Get the base file name (remove the directory path)
    base_file_name=$(basename "$file")
    # Remove the file extension
    file_name="${base_file_name%.*.*}"

    mkdir -p $folder_data/phylogenomics_analysis/trees/$set_name/seq_core/$file_name

    iqtree \
        -nt AUTO \
        -B 1000 \
         -s $folder_genomics/pangenome/$set_name/aligned_gene_sequences/$file_name.aln.fas \
        --prefix $folder_data/phylogenomics_analysis/trees/$set_name/seq_core/$file_name/$file_name
done

