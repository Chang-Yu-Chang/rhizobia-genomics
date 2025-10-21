#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script compares the ANI between the query and the reference fasta

# Prepare the list of query and reference
ref_list="$folder_genomics/taxonomy/list_refs.txt"
query_list="$folder_genomics/taxonomy/list_query.txt"

find "$folder_genomics/taxonomy/refs" -type f -name "*.fna" > "$ref_list"
find "$folder_genomics/fasta/genomes" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) > "$query_list"

outdir="$folder_genomics/taxonomy/fastani_results"
mkdir -p "$outdir"

# Split list_refs.txt into chunks of 100 lines each (numeric suffixes, zero-padded to 3 digits)
split -l 100 -d -a 3 "$ref_list" "${outdir}/ref_batch_" && \
for f in ${outdir}/ref_batch_*; do mv "$f" "$f.txt"; done

# Run FastANI for each query × batch combination
mamba activate fastani
i=0
while IFS= read -r query; do
  ((i++))
  qbase=$(basename "$query")
  qbase=${qbase%.*}  # safely strip any extension (.fna, .fasta, etc.)
  echo ""
  echo "=== [${i}] Processing query: $qbase ==="

  for ref_batch in ${outdir}/ref_batch_*.txt; do
    refbase=$(basename "$ref_batch" .txt)
    outfile="${outdir}/${qbase}_vs_${refbase}.txt"

    # Skip if output already exists
    if [[ -s "$outfile" ]]; then
      echo "Skipping $outfile (already exists)"
      continue
    fi

    echo ">>> Running: $qbase vs $refbase"
    start_time=$(date +%s)

    fastANI \
      --query "$query" \
      --refList "$ref_batch" \
      --threads 6 \
      --output "$outfile"

    end_time=$(date +%s)
    runtime=$(( (end_time - start_time) / 60 ))
    echo ">>> Done: $outfile (${runtime} min)"
  done
done < "$query_list"

echo ""
echo "✅ All FastANI comparisons finished successfully."
