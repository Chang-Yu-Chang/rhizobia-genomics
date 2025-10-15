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

# Install blast v2.14.1
# blast stands for Basic Local Alignment Search Tool, and can build a BLAST database with local sequences
mamba create -n blast
mamba activate blast
mamba install -y -c bioconda blast=2.14.1

# Install barrnap v0.9
# Barrnap is BAsic Rapid Ribosomal RNA Predictor.  predicts the location of ribosomal RNA genes in genomes
mamba create -n barrnap
mamba activate barrnap
mamba install -y -c bioconda barrnap=0.9
# Downalod the RefSeq 16S database
#https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/ # Proceed to Send to > Complete Record > File > FASTA > Sort by Default order
# Store it at ~/bioinformatics/16s/refseq_16s.fasta. Make it into a database using blast
mkdir -p ~/bioinformatics/16s/
mv ~/Downloads/sequence.fasta ~/bioinformatics/16s/refseq_16s.fasta
mamba activate blast
makeblastdb -in ~/bioinformatics/16s/refseq_16s.fasta -dbtype nucl

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
mamba install -y -c conda_forge biopython=1.84
mamba install -y -c conda_forge pandas
mamba install -y matplotlib

# Install gffread v0.12.7
mamba create -n gffread
mamba activate gffread
mamba install -y -c bioconda gffread=0.12.7

# Install defense finder v1.3.0
mamba create -n defense-finder
mamba activate defense-finder
mamba install -y -c bioconda defense-finder=1.3.0
# Install the model from macsyfinder
pip install -U mdmparis-defense-finder
defense-finder update

# Install padloc v2.0.0
mamba create -y -n padloc -c conda-forge -c bioconda -c padlocbio padloc=2.0.0
mamba activate padloc
# Download the latest database
padloc --db-update

# Install clonalframml
mamba create -n clonalframeml
mamba activate clonalframeml
mamba install -y -c bioconda clonalframeml=1.13

# Install networkx
mamba create -n networkx
mamba activate networkx
mamba install -y networkx

# Install carveme 1.6.4
# CarveMe: automated genome-scale metabolic model reconstruction
mamba create -n carveme
mamba activate carveme
mamba install -y -c bioconda carveme=1.6.4
mamba install -y -c bioconda diamond=2.1.12
mamba install -y -c conda-forge pyscipopt=5.5.0
mamba install -y pandas=2.2.3 # new pandas introduces syntax error

# Install cobra 0.21.0
mamba create -n cobra
mamba activate cobra
mamba install -y python=3.11
mamba install -y -c conda-forge cobra=0.29.1 numpy=1.23.5

# Install cogclassifier 2.0.0
# Classify prokaryote protein sequences into COG functional category.
mamba create -n cogclassifier
mamba activate cogclassifier
mamba install -y -c bioconda cogclassifier=2.0.0

# Install CheckM 1.2.4
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
