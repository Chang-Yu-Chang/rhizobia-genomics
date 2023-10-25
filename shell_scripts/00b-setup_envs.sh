# This script creates the mamba envs for various bioconda tools.
# In priniciple, one env per one tool

cd

# Install Nanocomp v1.23.1
# Nanocomp compares multiple runs of long read sequencing data and alignments.
mamba create -n nanocomp
mamba activate nanocomp
mamba install --yes -c bioconda nanocomp=1.23.1

# Install Nanoplot v.1.41.6
# Nanoplot is a plotting tool for long read sequencing data and alignments.
mamba create -n nanoplot
mamba activate nanoplot
mamba install --yes -c bioconda nanoplot=1.41.6

# Install Filtlong v0.2.1
# Filtlong is a tool for filtering long reads by quality.
mamba create -n filtlong
mamba activate filtlong
mamba install --yes -c bioconda filtlong=0.2.1

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

# # Install Bandage v0.8.1
# # Bandage (a Bioinformatics Application for Navigating De novo Assembly Graphs Easily)
# mamba create -n bandage
# mamba activate bandage
# mamba install bandage=0.8.1

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


# Install CheckM v1.2.2
# CheckM, an automated method for assessing the quality of a genome using a broader set of marker genes specific to the position of a genome within a reference genome tree and information about the collocation of these gene
mamba create -n checkm
mamba activate checkm
mamba install --yes -c bioconda checkm-genome=1.2.2

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
# Download a prefetch database for k=3.1
# # Download a Genbank LCA database for k=31
# mkdir -p ~/bioinformatics/sourmash
# cd ~/bioinformatics/sourmash/
# wget https://osf.io/4f8n3/download -O genbank-k31.lca.json.gz



