# Biogeographic and Genomic Signatures of Thermal Adaptation in Facultative Symbionts

**Chang-Yu Chang, Terrence Topping-Brown, Jazmine L. Rud, McCall B. Calvert, Gerardo Bencosme, and Corlett Wood**

This repository contains scripts and minimal data to reproduce all analyses and figures from the manuscript:

> **“Biogeographic and genomic signatures of thermal adaptation in facultative symbionts.”**


## Setup

This project uses **bioconda/miniconda** with **mamba** for bioinformatic tool management and **renv** for R package reproducibility.

### 1. Conda / Mamba

1. Install conda, miniforge, or mamba using the commands in `setup_conda.sh`.  
2. Once `mamba` is installed, run `setup_envs.sh` to create the individual bioconda environments.  
   - Each bioinformatic tool (e.g., *panaroo*, *mafft*, *fastANI*) resides in its own environment.

### 2. R Environment

The R scripts were executed under the following environment:

```
> sessionInfo()
R version 4.4.2 (2024-10-31)
Platform: x86_64-apple-darwin20
Running under: macOS Sequoia 15.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Taipei
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.4.2      BiocManager_1.30.25 tools_4.4.2    rstudioapi_0.17.1   renv_1.0.9   
```

To initialize `renv` and restore packages:

```r
> install.packages("renv")
> renv::restore()
```

# Data folders and files

1. `mapping/` has a mapping file
    - `isolates.csv` is the list of isolates/strains and their sequencing batches and ids
2. `genomics/`
    - `annotation` contains the Prokka annotation output, one folder per genome
    - `assembly` contains the assembly pipiline output, one folder per genome
    - `checkm` stores the checkm output
    - `fasta` contains the consolidated genome fasta that were uploaded to NCBI BioProject PRJNA1338806
    - `gff` has the consolidated genome gff
    - `ibd` stores the IBD analysis output
    - `pangenome` has the panaroo output, gene content analysis, and tree making
    - `qcs` has the consolidated assembly QC tables
    - `raw_reads` contains the fastq.gz files that were uploaded to NCBI BioProject PRJNA1338806
    - `taxonomy` has the taxonomy assignment results from the ANI-based approach
4. `phenotypes/` 
    - `sites/`: field sampling sites
    - `growth/`: growth traits

# Workflow

The shell scripts were executed on a 2021 iMac with Apple M1 chip 16GB memory and macOS version 14.7. 

## Reproducing figures and tables

The folder `plotting_scripts/` includes Rscripts each generates a figure or table in the manuscript.

To run all Rscripts, execute the following shell script:

```
for file in plotting_scripts/*.R; do
    Rscript -e "renv::activate('.'); source('$file')"
done
```

## Reproducing anaylsis from raw data

To reproduce analysis from raw data (raw reads and trait data), analysis, to the final figures and tables, execute the following shell script:

```
# --- Growth assays ---
Rscript -e "renv::activate('.'); source('phenotypes/growth/fit_gc.R')"

# --- Site climate data (DAYMET) ---
Rscript -e "renv::activate('.'); source('phenotypes/sites/extract_climate.R')"

# --- Assembly pipeline ---
cd genomics/assembly/
zsh assess_reads.sh
zsh denovo_assembly.sh
zsh assess_assemblies.sh
zsh manual_concat.sh
cd ../../
Rscript -e "renv::activate('.'); source('genomics/assembly/consolidate_qcs.R')"

# --- Annotation ---
cd genomics/annotation/
zsh annotate_genomes.sh

# --- Taxonomy ---
cd genomics/taxonomy/
zsh download_ncbi_genomes.sh
zsh ani.sh
cd ../../
Rscript -e "renv::activate('.'); source('genomics/taxonomy/consolidate_ani.R')"

# --- Pangenome ---
cd genomics/pangenome/
zsh pangenome.sh
cd ../../
Rscript -e "renv::activate('.'); source('genomics/pangenome/clean_gpa.R')"
cd genomics/pangenome/
zsh concatenate_alignment.sh
zsh compute_trees1.sh
cd ../../
Rscript -e "renv::activate('.'); 
source('genomics/pangenome/compute_trees2.R')
source('genomics/pangenome/curate_trees.R')
source('genomics/pangenome/rf_tree.R')"

# --- IBD ---
zsh blash.sh
Rscript -e "renv::activate('.'); 
source('genomics/ibd/clean_blash.R')
source('genomics/ibd/make_gpa_by_replicon.R')
source('genomics/ibd/compute_sccg_snps.R')
source('genomics/ibd/compute_gcv_diff.R')"
```
