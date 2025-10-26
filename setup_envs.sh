#!/usr/bin/env zsh

# This script creates the mamba envs for various bioconda tools
# In priniciple, one env per tool

cd
mkdir -p ~/bioinformatics

# Install filtlong v0.2.1
# filtlong is a tool for filtering long reads by quality.
mamba create -n filtlong
mamba activate filtlong
mamba install -y -c bioconda filtlong=0.2.1

# Install bioawk v1.0
# bioawk is BWK awk modified for biological data
mamba create -n bioawk
mamba activate bioawk
mamba install -y -c bioconda bioawk=1.0

# Install miniasm v0.3
# miniasm is a very fast OLC-based de novo assembler for noisy long reads
# This is used for drafting genomes from raw reads
mamba create -n miniasm
mamba activate miniasm
mamba install -y -c bioconda minimap2=2.26
mamba install -y -c bioconda miniasm=0.3
mamba install -y -c bioconda any2fasta=0.4.2

# Install flye assembler v2.9.2
# flye is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies
mamba create -n flye
mamba activate flye
mamba install -y -c bioconda flye=2.9.2

# Install medaka v1.8.0
# medaka is a tool to create consensus sequences and variant calls from nanopore sequencing data
mamba create -n medaka
mamba activate medaka
mamba install -y -c bioconda medaka=1.8.0

# Install quast v5.2.0
# quast evaluates genome/metagenome assemblies by computing various metrics.
mamba create -n quast
mamba activate quast
mamba install -y -c bioconda quast=5.2.0

# Install BUSCO v6.0.0
# BUSCO: Assessment of assembly completeness using Universal Single Copy Orthologs
mamba create -n busco
mamba activate busco
mamba install -y -c bioconda busco=6.0.0 python=3.10

# Install prokka v1.14.5
# prokka is for rapid prokaryotic genome annotation
mamba create -n prokka
mamba activate prokka
mamba install -y -c bioconda prokka=1.14.5

# Install ncbi-datasets v18.9.0
# ncbi-datasets is NCBI Datasets command-line tools
mamba create -y -n ncbi-datasets
mamba activate ncbi-datasets
mamba install -y -c bioconda ncbi-datasets-cli=18.9.0

# Install minimap2 v2.26
# minimap2 is a versatile pairwise aligner for genomic and spliced nucleotide sequences
mamba create -y -n minimap2
mamba activate minimap2
mamba install -y -c bioconda minimap2=2.26

# Install fastANI v1.31
# fastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
mamba create -n fastani
mamba activate fastani
mamba install -y -c bioconda fastani=1.31

# Install panaroo v1.3.4
# panaroo is A Bacterial Pangenome Analysis Pipeline that can call large structural variants
mamba create -n panaroo
mamba activate panaroo
mamba install -y -c bioconda panaroo=1.3.4

# Install iqtree v2.3.0
# iqtree is an Efficient phylogenomic software by maximum likelihood.
mamba create -n iqtree
mamba activate iqtree
mamba install -y -c bioconda iqtree=2.3.0

# Install biopython v1.84
mamba create -n biopython
mamba activate biopython
mamba install -y -c conda-forge biopython=1.84
mamba install -y -c conda-forge pandas=2.3.3
mamba install -y matplotlib

# Install gffread v0.12.7
mamba create -n gffread
mamba activate gffread
mamba install -y -c bioconda gffread=0.12.7

# Install checkM 1.2.4
mamba create -n checkm python=3.9
mamba activate checkm
mamba install -y -c bioconda numpy matplotlib pysam
mamba install -y -c bioconda hmmer prodigal pplacer
pip install checkm-genome
# Download reference data
curl -o ~/.checkm/checkm_data_2015_01_16.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
cd ~/.checkm
tar -xvzf checkm_data_2015_01_16.tar.gz
checkm data setRoot ~/.checkm

# Install bedtools v2.31.1
mamba create -n bedtools
mamba activate bedtools
mamba install -y -c bioconda bedtool=2.31.1

# Install mafft v7.525
mamba create -n mafft
mamba activate mafft
mamba install y -c bioconda mafft=7.525

# Install mafft v1.4.1
mamba create -n trimal
mamba activate trimal
mamba install y -c bioconda trimal=1.4.1
