#!/usr/bin/env zsh

# This script creates the mamba envs for various bioconda tools
# In priniciple, one env per tool

cd
mkdir -p ~/bioinformatics

# Install Filtlong v0.2.1
# Filtlong is a tool for filtering long reads by quality.
mamba create -n filtlong
mamba activate filtlong
mamba install --yes -c bioconda filtlong=0.2.1

# Install bioawk v1.0
mamba create -n bioawk
mamba activate bioawk
mamba install --yes -c bioconda bioawk=1.0

# Install miniasm v0.3
# Miniasm is a very fast OLC-based de novo assembler for noisy long reads
mamba create -n miniasm
mamba activate miniasm
mamba install --yes -c bioconda minimap2=2.26
mamba install --yes -c bioconda miniasm=0.3
mamba install --yes -c bioconda any2fasta=0.4.2

# Install Flye assembler v2.9.2
# Flye is a de novo assembler for single-molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies
mamba create -n flye
mamba activate flye
mamba install --yes -c bioconda flye=2.9.2

# Install Medaka v1.8.0
# medaka is a tool to create consensus sequences and variant calls from nanopore sequencing data. This task is performed using neural networks applied a pileup of individual sequencing reads against a draft assembly. It provides state-of-the-art results outperforming sequence-graph based methods and signal-based methods, whilst also being faster.
mamba create -n medaka
mamba activate medaka
mamba install --yes -c bioconda medaka=1.8.0



# Install ncbi-datasets v15.27.1
mamba create -y -n ncbi-datasets
mamba activate ncbi-datasets
mamba install --yes -c bioconda ncbi-datasets=15.27.1

# Install minimap2 v2.26
mamba create -y -n minimap2
mamba activate minimap2
mamba install --yes -c bioconda minimap2=2.26

# # Install BWA v0.7.17
# # BWA is a software package for mapping low-divergent sequences against a large reference genome,
# mamba create -y -n bwa
# mamba activate bwa
# mamba install --yes -c bioconda bwa=0.7.17

# Install Samtools v1.18
# samtools are Tools for dealing with SAM, BAM and CRAM files
mamba create -y -n samtools
mamba activate samtools
mamba install --yes -c bioconda samtools=1.18

# Install Snippy v4.6.0
# Snippy is Rapid haploid variant calling and core genome alignment
mamba create -y -n snippy
mamba activate snippy
mamba install --yes -c bioconda snippy=4.6.0
mamba install --yes -c bioconda vcflib=1.0.1 # Because vcflib 1.0.2 breaks snippy. Downgrade it according to https://github.com/tseemann/snippy/issues/561

# # Install bcftools 1.18
# # BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF
# mamba create -y -n bcftools
# mamba activate bcftools
# mamba install --yes -c bioconda bcftools=1.18
#
# # Install vcftools 0.1.16
# # vcftools is A set of tools written in Perl and C++ for working with VCF files.
# mamba create -y -n vcftools
# mamba activate vcftools
# mamba install --yes -c bioconda vcftools=0.1.16

# Install sniffles2 v2.2
# A fast structural variant caller for long-read sequencing, Sniffles2 accurately detect SVs on germline, somatic and population-level for PacBio and Oxford Nanopore read data.
mamba create -y -n sniffles
mamba activate sniffles
mamba install --yes -c bioconda sniffles=2.2


# Install Mash v2.3
# Mash is a tool commonly used for fast and memory-efficient sequence similarity estimation and taxonomy classification.
mamba create -n mash
mamba activate mash
mamba install --yes -c bioconda mash=2.3
# Download the RefSeq database. It includes both genomes and plasmids
mkdir -p ~/bioinformatics/mash
cd ~/bioinformatics/mash/
wget https://gembox.cbcb.umd.edu/mash/refseq.genomes%2Bplasmid.k21s1000.msh -O refseq.genomes%2Bplasmid.k21.s1000.msh

# Install Sourmash v4.8.4
# Quickly search, compare, and analyze genomic and metagenomic data sets
mamba create -n sourmash
mamba activate sourmash
mamba install --yes -c bioconda sourmash=4.6.1
# Download a prefetch GTDB R8 genomic database for k=31
mkdir -p ~/bioinformatics/sourmash
cd ~/bioinformatics/sourmash/
# Install zip that is used for sourmash gather
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.zip -O gtdb-rs214-k31.zip


# # genomic representatives; this is 4.4GB with 85,205 species-level genomes. This took ~7min
# wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.sbt.zip -O gtdb-rs214-reps.k31.sbt.zip
# # genomic all genomes; this is 23GB with 402,709 genomes. This took 40min
# wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.sbt.zip -O gtdb-rs214-k31.sbt.zip
# # download species-level lineages
# wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214.lineages.csv.gz -O gtdb-rs214.lineages.csv.gz
# gunzip gtdb-rs214.lineages.csv.gz
#
#wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2022.03/genbank-2022.03-bacteria-k31.zip -O genbank-2022.03-bacteria-k31.zip
# Create a Sequence Bloom Tree (SBT) database from the download database
#sourmash index -k 31 genbank-2022.03-bacteria-k31.sbt.zip genbank-2022.03-bacteria-k31.zip
# `sourmash index` to create a Sequence Bloom Tree (SBT) that can be quickly searched on disk; this is the same format in which we provide GenBank and other databases
# This will take a while
#wget https://osf.io/4f8n3/download -O genbank-k31.lca.json.gz

# Install barrnap v0.9
# Barrnap is BAsic Rapid Ribosomal RNA Predictor.  predicts the location of ribosomal RNA genes in genomes
mamba create -n barrnap
mamba activate barrnap
mamba install --yes -c bioconda barrnap=0.9
# Downalod the RefSeq 16S database
#https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/ # Proceed to Send to > Complete Record > File > FASTA > Sort by Default order
# Make it into a database
#makeblastdb -in $refseq_16s_db -dbtype nucl

# Install blast v2.14.1
mamba create -n blast
mamba activate blast
mamba install --yes -c bioconda blast=2.14.1


# Install prokka v1.14.5
mamba create -n prokka
mamba activate prokka
mamba install --yes -c bioconda prokka=1.14.5











# Install Bakta v1.6.1
# bakta, a new command-line software tool for the robust, taxon-independent, thorough and, nonetheless, fast annotation of bacterial genomes.
mamba create -n bakta
mamba activate bakta
mamba install --yes -c bioconda bakta=1.6.1
# Download the mandatory database
bakta_db list # show available database
mkdir -p ~/bioinformatics/bakta
bakta_db download --output ~/bioinformatics/bakta/
# This database is on zenodo https://zenodo.org/records/7669534 and it's ~50 GB

# Install Quast v5.2.0
# QUAST stands for QUality ASsessment Tool. It evaluates genome/metagenome assemblies by computing various metrics. The current QUAST toolkit includes the general QUAST tool for genome assemblies, MetaQUAST, the extension for metagenomic datasets, QUAST-LG, the extension for large genomes (e.g., mammalians), and Icarus, the interactive visualizer for these tools.
mamba create -n quast
mamba activate quast
mamba install --yes -c bioconda quast=5.2.0
# The default QUAST package does not include:
# * GRIDSS (needed for structural variants detection)
# * SILVA 16S rRNA database (needed for reference genome detection in metagenomic datasets)
# * BUSCO tools and databases (needed for searching BUSCO genes) -- works in Linux only!
# To be able to use those, please run
quast-download-gridss
quast-download-silva
quast-download-busco

# Install Busco v5.4.7
# BUSCO provides measures for quantitative assessment of genome assembly, gene set, and transcriptome completeness based on evolutionarily informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB.
mamba create -n busco
mamba activate busco
mamba install --yes -c bioconda busco=5.4.7



# Install FastANI v1.31
# FastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
mamba create -n fastani
mamba activate fastani
mamba install --yes -c bioconda fastani=1.31

# Install Roary v3.13.0
# Roary is a tool designed to quickly build large-scale pan genomes for prokaryote populations
mamba create -n roary
mamba activate roary
mamba install --yes -c bioconda roary=3.13.0
roary -h
# roary may plot using R and ggplot2
mamba install --yes -c conda-forge r-base=4.2.0
mamba install --yes -c conda-forge r-ggplot2=3.4.4
# for plotting tree
mamba install --yes -c bioconda fasttree=2.1.11
# an additional plotting script using the gene presence and absence table
wget https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py -O ~/bioinformatics/roary/roary_plots.py
# this requires the following python packages in the env
mamba install --yes -c conda-forge matplotlib=3.8.1
mamba install --yes -c conda-forge seaborn=0.13.0
mamba install --yes -c conda-forge biopython=1.81

# Install anvi'o v8
# Anvi'o is an open-source, community-driven analysis and visualization platform for microbial 'omics.
# Set up environment
mamba create -y -n anvio-8 python=3.10
mamba activate anvio-8
mamba install -y -c conda-forge -c bioconda python=3.10 \
        sqlite prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript
#mamba install -y -c bioconda fastani
# Download the python source code package
mkdir -p ~/bioinformatics/anvio
curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz --output ~/bioinformatics/anvio/anvio-8.tar.gz
# For macOS that requires more up to date c-compiler
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
# pip install
cd ~/bioinformatics/anvio
pip3 install anvio-8.tar.gz # make sure the pip3 is under the anvio environment
# Check installation. If everything goes smoothly, your browser should pop-up and show you an anvi’o interactive interface
anvi-self-test --suite mini
# Setup key resources
anvi-setup-scg-taxonomy # to setup SCG taxonomy data using GTDB genomes.
anvi-setup-ncbi-cogs, # to setup NCBI’s COG database for quick annotation of genes with functions,
anvi-setup-kegg-data, # so anvi-estimate-metabolism and anvi-reaction-network find the database of KEGG orthologs ready when you need it.
anvi-self-test --suite pangenomics # to see if everything is order, especially if you plan to use anvi’o for pangenomics.
# Install package that converts prokka annotation to anvio compatible format
cd ~/bioinformatics/anvio
wget https://raw.githubusercontent.com/karkman/gff_parser/master/gff_parser.py -O gff_parser.py
pip3 install gffutils

# Install PlasmidFinder
mamba create -y -n plasmidfinder
mamba activate plasmidfinder
# Go to wanted location for plasmidfinder
cd ~/bioinformatics
mamba install --yes -c bioconda plasmidfinder=2.1.6
# Build container
docker build -t plasmidfinder .
docker run --rm -it --entrypoint=/test/test.sh plasmidfinder
# Install PlasmidFinder database with executable kma_index program
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
cd plasmidfinder_db
rm -rf .git # remove the git
PLASMID_DB=$(pwd)
python3 INSTALL.py kma_index


# Install Freebayes v1.3.6
mamba create -y -n freebayes
mamba activate freebayes
mamba install --yes -c bioconda freebayes=1.3.6


# Install ggcaller v1.3.2
# ggCaller traverses Bifrost graphs constructed from bacterial genomes to identify putative gene sequences, known as open reading frames (ORFs).
# ggCaller incorporates Balrog to filter ORFs to improve specificity of calls and Panaroo for pangenome analysis and quality control.
mamba create -y -n ggcaller
mamba activate ggcaller
mamba install --yes -c bioconda ggcaller=1.3.2

# Install plink2 v2.00a5
mamba create -y -n plink2
mamba activate plink2
mamba install --yes -c bioconda plink2=2.00a5






# # Install Nanocomp v1.23.1
# # Nanocomp compares multiple runs of long read sequencing data and alignments.
# mamba create -n nanocomp
# mamba activate nanocomp
# mamba install --yes -c bioconda nanocomp=1.23.1

# # Install Nanoplot v.1.41.6
# # Nanoplot is a plotting tool for long read sequencing data and alignments.
# mamba create -n nanoplot
# mamba activate nanoplot
# mamba install --yes -c bioconda nanoplot=1.41.6

# # Install Bandage v0.8.1
# # Bandage (a Bioinformatics Application for Navigating De novo Assembly Graphs Easily)
# mamba create -n bandage
# mamba activate bandage
# mamba install bandage=0.8.1

# # Install CheckM v1.2.2
# # CheckM, an automated method for assessing the quality of a genome using a broader set of marker genes specific to the position of a genome within a reference genome tree and information about the collocation of these gene
# mamba create -n checkm
# mamba activate checkm
# mamba install --yes -c bioconda checkm-genome=1.2.2


# # Install Scoary v1.6.16
# # Scoary is a tool that uses Roary's gene_presence_absence.csv file and a user-created traits file to calculate associations between accessory genome genes and traits. It generates a list of genes ranked by the strength of association for each trait.
# mamba create -n scoary
# mamba activate scoary
# mamba install --yes -c bioconda scoary=1.6.16
# scoary -h
#
# # Install Artemis  v18.2.0
# # Artemis is a free genome browser and annotation tool that allows visualisation of sequence features, next generation data and the results of analyses within the context of the sequence, and also its six-frame translation.
# mamba create -n artemis
# mamba activate artemis
# mamba install --yes -c bioconda artemis=18.2.0


# # Clone and enter the plasmidfinder directory
# git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
# rm -rf plasmidfinder/.git # remove the git
# mamba env update -f conda/meta.yaml

#Please run download-db.sh to download the PlasmidFinder database to /Users/cychang/miniconda3/envs/intel_env/envs/plasmidfinder/share/plasmidfinder-2.1.6/database.

# Install plaSquid
# # plaSquid is a Nextflow pipeline for plasmid detection and classification in genomic and metagenomic data.
# mkdir ~/bioinformatics/plasquid
# cd ~/bioinformatics/plasquid
# git clone --depth=1 https://github.com/mgimenez720/plaSquid/
# rm -rf plaSquid/.git # remove the git
# # Create mamba env
# cd ~/bioinformatics/plasquid/plaSquid
# mamba create -y -n plasquid
# mamba activate plasquid
# mamba env update -f environments/plaSquid.yml
# mamba search -c conda-forge procpsng procps-ng
#
# # Install SDK
# curl -s "https://get.sdkman.io" | bash
# # Install java
# sdk install java 17.0.6-amzn
# # Install Nextflow
# take ~/bioinformatics/nextflow
# wget -qO- https://get.nextflow.io | bash
# # Test
# cd ~/bioinformatics/plasquid/plaSquid
# chmod 755 ~/bioinformatics/nextflow/nextflow
# ~/bioinformatics/nextflow/nextflow run main.nf --contigs testdata/test.fasta --outdir testdata

# # Install PlasForest
# # PlasForest is A random forest classifier of contigs to identify contigs of plasmid origin in contig and scaffold genomes.
# # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04270-w
# mamba create -y -n plasforest python=3.8 # scikit-learn==0.22.2 requires python <= 3.8
# mamba activate plasforest
# mamba install --yes -c bioconda biopython=1.79
# mamba install --yes -c bioconda numpy=1.21
# mamba install --yes -c bioconda pandas=2.0.1
# mamba install --yes -c bioconda joblib=1.3.0
# mamba install --yes -c bioconda scikit-learn=0.22.2.post1
# mamba install --yes -c bioconda blast=2.10.1
# # Download the code
# mkdir ~/bioinformatics/plasforest
# cd ~/bioinformatics/plasforest
# git clone --depth=1 https://github.com/leaemiliepradier/PlasForest
# rm -rf PlasForest/.git # remove the git
# # Untar the random forest classifier:
# cd PlasForest/
# tar -zxvf plasforest.sav.tar.gz
# # Download a database of plasmid sequences (2.5GB)
# chmod 755 database_downloader.sh
# ./database_downloader.sh
# # Testing the installation. To test that PlasForest has been correctly installed, you can run the following script:
# chmod 755 ./test_plasforest.sh
# ./test_plasforest.sh

















