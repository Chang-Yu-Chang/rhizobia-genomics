#!/usr/bin/env zshs
source ~/.zshrc
source 00-env_vars.sh

# This script assigns taxonomy to the de novo assembled genomes and contigs

cd $folder_shell
refseq_db="/Users/cychang/bioinformatics/mash/refseq.genomes+plasmid.k21.s1000.msh"
gtdb_db="/Users/cychang/bioinformatics/sourmash/gtdb-rs214-k31.zip"
refseq_16s_db="/Users/cychang/bioinformatics/16s/refseq_16s.fasta"
blast_genomes_db="$folder_genomics/blast_db/genomes"

# Make customized genome database
zsh 04a-make_db.sh $folder_genomics

mkdir -p $folder_genomics/taxonomy

for i in {1..38}
do
    echo $genome_ids[$i]
    genome_fa="$folder_genomics/genomes/$genome_ids[$i].fasta"
    mkdir -p $folder_genomics/taxonomy/$genome_ids[$i]

    # Compare strains to database
    mkdir -p $folder_genomics/taxonomy/$genome_ids[$i]/sourmash
    zsh 04c-sourmash.sh \
        $genome_fa \
        $folder_genomics/taxonomy/$genome_ids[$i]/sourmash \
        $gtdb_db
    
    # Compare the k-mer signature among genomes
    mkdir -p $folder_genomics/taxonomy/$genome_ids[$i]/kmer
    zsh 04f-genome_kmer.sh \
        $genome_fa \
        $folder_genomics/taxonomy/$genome_ids[$i]/kmer/genome.sig




#     # Estimate genome distance via mash
#     mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/mash"
#     zsh 04b-mash.sh \
#         $genome_fa \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/mash" \
#         $refseq_db
#
#     # Compare genomes via sourmash
#     mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/sourmash"
#     zsh 04c-sourmash.sh \
#         $genome_fa \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/sourmash" \
#         $gtdb_db
#
#     # Extract 16S rRNA from genome and blast
#     mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s"
#     zsh 04d-blast_16s.sh \
#         $genome_fa \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/rrna.fasta" \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/rrna.txt" \
#         $refseq_16s_db \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/blast.txt" \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/16s/taxonomy.txt"
#
#     # Blast genomes to a customized database of meliloti and medicae
#     mkdir -p "$folder_genomes/$sample_ids[$i]/04-taxonomy/blast_genome"
#     zsh 04e-blast_genome.sh \
#         $genome_fa \
#         $blast_genomes_db \
#         "$folder_genomes/$sample_ids[$i]/04-taxonomy/blast_genome/taxonomy.txt"

done























