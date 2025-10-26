#!/usr/bin/env zsh
# Extract recA and dnaE2 from annotated genomes,
# check for duplicates, and build phylogenies

source ~/.zshrc
source ../env_vars.sh

gff_dir=$folder_genomics/gff
fna_dir=$folder_genomics/fasta/genomes
out_dir=$folder_genomics/pangenome/recA_dnaE2
threads=8

mkdir -p $out_dir
mkdir -p $out_dir/fastas
cp ${fna_dir}/*.fasta $out_dir/fastas/

# Extracting recA and dnaE2 CDSs
mamba activate bedtools
mkdir -p "$out_dir"/{fastas,beds,seqs}
cd $out_dir/beds/

for gff in ${gff_dir}/*.gff; do
    base=$(basename "$gff" .gff)

    # Support either .fna or .fasta input
    if [[ -f "$out_dir/fastas/${base}.fna" ]]; then
        fasta="$out_dir/fastas/${base}.fna"
    elif [[ -f "$out_dir/fastas/${base}.fasta" ]]; then
        fasta="$out_dir/fastas/${base}.fasta"
    else
        echo "[warning] FASTA not found for $base — skipping."
        continue
    fi

    for gene in recA dnaE2; do
        awk -v g="$gene" '$3=="CDS" && $0~g {print $0}' "$gff" > "${base}_${gene}.gff" || true
        if [[ -s "${base}_${gene}.gff" ]]; then
            awk -v OFS="\t" '{print $1,$4-1,$5,$9,$6,$7}' "${base}_${gene}.gff" > "${base}_${gene}.bed"
            bedtools getfasta -fi "$fasta" -bed "${base}_${gene}.bed" -s -name > "$out_dir/seqs/${base}_${gene}.fa"
        else
            echo "[info] No $gene hits in $base"
        fi
    done
done

# Concatenating all sequences
cd $out_dir/seqs
cat *_recA.fa > all_recA.fa
cat *_dnaE2.fa > all_dnaE2.fa
sed -E 's/ .*//g' -i all_recA.fa all_dnaE2.fa

# Creating BLAST databases
mamba activate blast
cd $out_dir/seqs
makeblastdb -in all_recA.fa -dbtype nucl -out recA_db
makeblastdb -in all_dnaE2.fa -dbtype nucl -out dnaE2_db

# Running self-BLAST for within-genome identity
blastn -query all_recA.fa -db recA_db -out $out_dir/recA_self.tsv \
  -outfmt "6 qseqid sseqid pident length qlen slen" -num_threads $threads
blastn -query all_dnaE2.fa -db dnaE2_db -out $out_dir/dnaE2_self.tsv \
  -outfmt "6 qseqid sseqid pident length qlen slen" -num_threads $threads

# Multiple alignment and tree building
mamba activate mafft
cd $out_dir/seqs
mafft --auto --thread $threads all_recA.fa > all_recA_aln.fa
mafft --auto --thread $threads all_dnaE2.fa > all_dnaE2_aln.fa

mamba activate trimal
trimal -in all_recA_aln.fa -out all_recA_trim.fa -automated1
trimal -in all_dnaE2_aln.fa -out all_dnaE2_trim.fa -automated1

mamba activate iqtree
iqtree -s all_recA_trim.fa -m MFP -bb 1000 -nt $threads -pre $out_dir/seqs/recA_tree
iqtree -s all_dnaE2_trim.fa -m MFP -bb 1000 -nt $threads -pre $out_dir/seqs/dnaE2_tree
