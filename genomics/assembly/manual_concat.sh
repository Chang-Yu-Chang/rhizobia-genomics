#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script manual concatenates the contigs based on blast result
# g20: concatenate contigs 3, 7, 9
# g24: concatenate contigs 10, 17, 9, 14, 18

mamba activate seqtk

# g20
dir=$folder_genomics/assembly/g20
echo ">contig_3" > $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_3\ncontig_7\ncontig_9") | grep -v "^>" >> $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_1\ncontig_2") >> $dir/final_genome.fasta
# ensure the output has the same line length
awk '/^>/ {print; getline; i=0; while (i < length($0)) { print substr($0, i+1, 60); i+=60; } }' $dir/final_genome.fasta > $dir/final_genome_fixed.fasta
cat $dir/final_genome_fixed.fasta | grep '>'

cp $dir/final_genome_fixed.fasta $folder_genomics/fasta/genomes/g20.fasta

# g24
dir=$folder_genomics/assembly/g24
echo ">contig_10" > $dir/final_genome.fasta # chromosome
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_10\ncontig_17\ncontig_9\ncontig_14\ncontig_18") | grep -v "^>" >> $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_1\ncontig_3") >> $dir/final_genome.fasta
awk '/^>/ {print; getline; i=0; while (i < length($0)) { print substr($0, i+1, 60); i+=60; } }' $dir/final_genome.fasta > $dir/final_genome_fixed.fasta
cat $dir/final_genome.fasta | grep '>'

cp $dir/final_genome_fixed.fasta $folder_genomics/fasta/genomes/g24.fasta
