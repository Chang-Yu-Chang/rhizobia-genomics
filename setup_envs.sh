#!/usr/bin/env zsh

# This script creates the mamba envs for various bioconda tools
# In priniciple, one env per tool

cd
mkdir -p ~/bioinformatics

# Install filtlong v0.2.1
# filtlong is a tool for filtering long reads by quality.
mamba create -n filtlong
mamba activate filtlong
mamba install --yes -c bioconda filtlong=0.2.1

# Install bioawk v1.0
# bioawk is BWK awk modified for biological data
mamba create -n bioawk
mamba activate bioawk
mamba install --yes -c bioconda bioawk=1.0

# Install miniasm v0.3
# miniasm is a very fast OLC-based de novo assembler for noisy long reads
# This is used for drafting genomes from raw reads
mamba create -n miniasm
mamba activate miniasm
mamba install --yes -c bioconda minimap2=2.26
mamba install --yes -c bioconda miniasm=0.3
mamba install --yes -c bioconda any2fasta=0.4.2

# Install flye assembler v2.9.2
# flye is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies
mamba create -n flye
mamba activate flye
mamba install --yes -c bioconda flye=2.9.2

# Install medaka v1.8.0
# medaka is a tool to create consensus sequences and variant calls from nanopore sequencing data
mamba create -n medaka
mamba activate medaka
mamba install --yes -c bioconda medaka=1.8.0

# Install quast v5.2.0
# quast evaluates genome/metagenome assemblies by computing various metrics.
mamba create -n quast
mamba activate quast
mamba install --yes -c bioconda quast=5.2.0

# Install checkm v1.2.2
# checkm assess the quality of microbial genomes recovered from isolates, single cells, and metagenomes.
mamba create -n checkm --yes python=3.9
mamba activate checkm
mamba install --yes -c bioconda numpy matplotlib pysam hmmer prodigal
pip3 install checkm-genome==1.2.2
# Download checkm reference data at https://data.ace.uq.edu.au/public/CheckM_databases
# Uncompress the download zip
mv ~/Downloads/checkm_data_2015_01_16 ~/bioinformatics/checkm
checkm data setRoot ~/bioinformatics/checkm/checkm_data_2015_01_16
# Install pplacer binary files from https://github.com/matsen/pplacer/releases/tag/v1.1.alpha17
# Move the pplacer bin to system path
mv Downloads/pplacer-Darwin-v1.1.alpha17-6-g5cecf99/* ~/miniconda3/envs/intel_env/envs/checkm/bin/
# Check installation
#checkm test ~/checkm_test_result

# Install prokka v1.14.5
# prokka is for rapid prokaryotic genome annotation
mamba create -n prokka
mamba activate prokka
mamba install --yes -c bioconda prokka=1.14.5

# Install ncbi-datasets v15.27.1
# ncbi-datasets is NCBI Datasets command-line tools
mamba create -y -n ncbi-datasets
mamba activate ncbi-datasets
mamba install --yes -c bioconda ncbi-datasets=15.27.1

# Install blast v2.14.1
# blast stands for Basic Local Alignment Search Tool, and can build a BLAST database with local sequences
mamba create -n blast
mamba activate blast
mamba install --yes -c bioconda blast=2.14.1

# Install sourmash v4.8.4
# Sourmash is a quickly search, compare, and analyze genomic and metagenomic data sets
mamba create -n sourmash
mamba activate sourmash
mamba install --yes -c bioconda sourmash=4.6.1
# Download a prefetch GTDB R8 genomic database for k=31
mkdir -p ~/bioinformatics/sourmash
cd ~/bioinformatics/sourmash/
# Install zip that is used for sourmash gather
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.zip -O gtdb-rs214-k31.zip

# Install barrnap v0.9
# Barrnap is BAsic Rapid Ribosomal RNA Predictor.  predicts the location of ribosomal RNA genes in genomes
mamba create -n barrnap
mamba activate barrnap
mamba install --yes -c bioconda barrnap=0.9
# Downalod the RefSeq 16S database
#https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/ # Proceed to Send to > Complete Record > File > FASTA > Sort by Default order
# Make it into a database
#makeblastdb -in $refseq_16s_db -dbtype nucl

# Install minimap2 v2.26
# minimap2 is a versatile pairwise aligner for genomic and spliced nucleotide sequences
mamba create -y -n minimap2
mamba activate minimap2
mamba install --yes -c bioconda minimap2=2.26

# Install vcftools 0.1.16
# vcftools is A set of tools written in Perl and C++ for working with VCF files.
mamba create -y -n vcftools
mamba activate vcftools
mamba install --yes -c bioconda vcftools=0.1.16

# Install snippy v4.6.0
# snippy is Rapid haploid variant calling and core genome alignment
mamba create -y -n snippy
mamba activate snippy
mamba install --yes -c bioconda snippy=4.6.0
mamba install --yes -c bioconda vcflib=1.0.1 # Because vcflib 1.0.2 breaks snippy. Downgrade it according to https://github.com/tseemann/snippy/issues/561

# Install fastANI v1.31
# fastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
mamba create -n fastani
mamba activate fastani
mamba install --yes -c bioconda fastani=1.31

# Install panaroo v1.3.4
# panaroo is A Bacterial Pangenome Analysis Pipeline that can call large structural variants
mamba create -n panaroo
mamba activate panaroo
mamba install --yes -c bioconda panaroo=1.3.4

# Install iqtree v2.3.0
# iqtree is an Efficient phylogenomic software by maximum likelihood.
mamba create -n iqtree
mamba activate iqtree
mamba install --yes -c bioconda iqtree=2.3.0

















# Install mash v2.3
# Mash is a tool commonly used for fast and memory-efficient sequence similarity estimation and taxonomy classification.
mamba create -n mash
mamba activate mash
mamba install --yes -c bioconda mash=2.3
# Download the RefSeq database. It includes both genomes and plasmids
mkdir -p ~/bioinformatics/mash
cd ~/bioinformatics/mash/
wget https://gembox.cbcb.umd.edu/mash/refseq.genomes%2Bplasmid.k21s1000.msh -O refseq.genomes%2Bplasmid.k21.s1000.msh



# Install samtools v1.18
# samtools are Tools for dealing with SAM, BAM and CRAM files
mamba create -y -n samtools
mamba activate samtools
mamba install --yes -c bioconda samtools=1.18

# Install bcftools 1.18
# BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF
mamba create -y -n bcftools
mamba activate bcftools
mamba install --yes -c bioconda bcftools=1.18


# Install sniffles2 v2.2
# sniffles is a fast structural variant caller for long-read sequencing, Sniffles2 accurately detect SVs on germline, somatic and population-level for PacBio and Oxford Nanopore read data.
mamba create -y -n sniffles
mamba activate sniffles
mamba install --yes -c bioconda sniffles=2.2


# Install roary v3.13.0
# roary is a tool designed to quickly build large-scale pan genomes for prokaryote populations
mamba create -n roary
mamba activate roary
mamba install --yes -c bioconda roary=3.13.0



# Install seqtk v1.4
# seqtk is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format.
mamba create -n seqtk
mamba activate seqtk
mamba install --yes -c bioconda seqtk=1.4



# Install pyseer v1.3.11
# Sequence Element Enrichment Analysis (SEER)
# pyseer uses linear models with fixed or mixed effects to estimate the effect of genetic variation in a bacterial population on a phenotype of interest, while accounting for potentially very strong confounding population structure. This allows for genome-wide association studies (GWAS) to be performed in clonal organisms such as bacteria and viruses.
mamba create -n pyseer
mamba activate pyseer
mamba install --yes -c bioconda pyseer=1.3.11

# Install bwa v0.7.17
# BWA is a software package for mapping low-divergent sequences against a large reference genome,
mamba create -y -n bwa
mamba activate bwa
mamba install --yes -c bioconda bwa=0.7.17

# Install DAGchainer vr120920-3
# DAGchainer identifies syntenic regions.
mamba create -y -n dagchainer
mamba activate dagchainer
mamba install --yes -c bioconda dagchainer=r120920-3
