#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script downloads the Rhizobiales Order type strain genomes from NCBI

mamba activate ncbi-datasets

# Download genome metadata
datasets summary genome taxon "Rhizobiales" \
    --reference --as-json-lines | dataformat tsv genome \
    --fields organism-name,accession,assminfo-name,annotinfo-name,annotinfo-release-date \
    > $folder_genomics/taxonomy/rhizobiales.txt

# Download genomes
folder_rhi=$folder_genomics/taxonomy/rhizobiales
datasets download genome taxon "Rhizobiales" --reference --include genome --dehydrated --filename $folder_rhi.zip
unzip -o $folder_rhi.zip -d $folder_rhi
datasets rehydrate --directory $folder_rhi

# Collect all FASTA files into one folder
folder_refs=$folder_genomics/taxonomy/refs
mkdir -p $folder_refs
find $folder_rhi/ncbi_dataset/data -type f -name "*.fna" -exec cp {} $folder_refs \;

find "$folder_rhi/ncbi_dataset/data" -type f -name "*.fna" | while read -r file; do
  filename=$(basename "$file")
  dest="$folder_refs/$filename"

  if [[ -f "$dest" ]]; then
    echo "Skipping $filename (already exists in $folder_refs)"
  else
    echo "Copying $filename → $folder_refs"
    cp "$file" "$dest"
  fi
done


# datasets summary genome taxon "Sinorhizobium" \
#     --reference --as-json-lines | dataformat tsv genome \
#     --fields organism-name,accession,assminfo-name,annotinfo-name,annotinfo-release-date \
#     > $folder_genomics/taxonomy/type_strains.txt
# datasets summary genome taxon "Ensifer" \
#     --reference --as-json-lines | dataformat tsv genome \
#     --fields organism-name,accession,assminfo-name,annotinfo-name,annotinfo-release-date \
#     | tail -n +2 >> $folder_genomics/taxonomy/type_strains.txt

# folder_sino=$folder_genomics/taxonomy/sinorhizobium
# datasets download genome taxon "Sinorhizobium" --reference --include genome --dehydrated --filename $folder_sino.zip
# unzip -o $folder_sino.zip -d $folder_sino
# datasets rehydrate --directory $folder_sino
#
# folder_ensi=$folder_genomics/taxonomy/ensifer
# datasets download genome taxon "Ensifer" --reference --include genome --dehydrated --filename $folder_ensi.zip
# unzip -o $folder_ensi.zip -d $folder_ensi
# datasets rehydrate --directory $folder_ensi
#
# find $folder_sino/ncbi_dataset/data -type f -name "*.fna" -exec cp {} $folder_type \;
# find $folder_ensi/ncbi_dataset/data -type f -name "*.fna" -exec cp {} $folder_type \;

# Clean
# rm $folder_genomics/taxonomy/sinorhizobium.zip
# rm $folder_genomics/taxonomy/ensifer.zip
# rm -r $folder_genomics/taxonomy/sinorhizobium
# rm -r $folder_genomics/taxonomy/ensifer
