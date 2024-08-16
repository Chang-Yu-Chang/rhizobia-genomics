#!/usr/bin/env zsh
source ~/.zshrc
source ../env_vars.sh

# This script manual concatenates the contigs based on blast result
# need to concatenate contigs: g20, g24
# need to concatenate contigs and remove small
# need to divide contigs: g42

mamba activate seqtk

# g20
dir=$folder_genomics/assembly/g20
echo ">contig_3" > $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_3\ncontig_7\ncontig_9") | grep -v "^>" >> $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_1\ncontig_2") >> $dir/final_genome.fasta
cat $dir/final_genome.fasta | grep '>'

cp $dir/final_genome.fasta $folder_genomics/fasta/genomes/g20.fasta

# g24
dir=$folder_genomics/assembly/g24
echo ">contig_10" > $dir/final_genome.fasta # chromosome
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_10\ncontig_17\ncontig_9\ncontig_14\ncontig_18") | grep -v "^>" >> $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_1\ncontig_3") >> $dir/final_genome.fasta
cat $dir/final_genome.fasta | grep '>'

cp $dir/final_genome.fasta $folder_genomics/fasta/genomes/g24.fasta

# g28
dir=$folder_genomics/assembly/g28
echo ">contig_12" > $dir/final_genome.fasta # chromosome
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_12\ncontig_17\ncontig_35\ncontig_79\ncontig_27\ncontig_78") | grep -v "^>" >> $dir/final_genome.fasta
echo ">contig_5" >> $dir/final_genome.fasta # psymA
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_5\ncontig_2\ncontig_13\ncontig_3\ncontig_9")  | grep -v "^>" >> $dir/final_genome.fasta
seqtk subseq $dir/medaka/consensus.fasta <(echo -e "contig_14") >> $dir/final_genome.fasta # psymB
cat $dir/final_genome.fasta | grep '>'

cp $dir/final_genome.fasta $folder_genomics/fasta/genomes/g28.fasta
