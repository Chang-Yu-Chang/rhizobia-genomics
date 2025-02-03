# General

Scripts and minimal data for reproducing the anaylsis, tables, and figures presented in the manuscript

Chang, Chang-Yu, Terrence Topping-Brown, Jazmine L. Rud, Mccall B. Calvert, Gerardo Bencosme, Linda J. Robinson, and Corlett Wood. 2025. “Testing for Divergence in a Plant Symbiont across Two Natural Environmental Gradients.” bioRxiv. https://doi.org/10.1101/2025.01.14.632453.

# Setup

This repository relies heavily on `conda`/`mamba` for managing bioinformatic tools and on `renv` for R packages.

### conda/mamba

- Follow commands in `setup_conda.sh` to install conda/miniforge/mamba on your terminal.

- Once you have installed `mamba`, follow commands in `setup_envs.sh` to install bioconda packages. In general, each bioconda tool has one mamba virtual environment. This includes the following packages at version:

filtlong v0.2.1
bioawk v1.0
miniasm v0.3
flye v2.9.2
medaka v1.8.0
quast v5.2.0
busco v5.7.1
prokka v1.14.5
ncbi-datasets v15.27.1
blast v2.14.1
sourmash v4.8.4 and a prefetch GTDB R8 genomic database for k=31
barrnap v0.9 and the RefSeq 16S database
minimap2 v2.26
fastani v1.31
panaroo v1.3.4
iqtree v2.3.0
biopython v1.84
gffread v0.12.7

### renv

The Rscripts are executed under the following R environment

```
> sessionInfo()
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
[1] BiocManager_1.30.25 compiler_4.3.2      tools_4.3.2         rstudioapi_0.17.1  
[5] renv_1.0.11       
```

For first time `renv` user, install renv and restore the packages recorded in `renv.lock`

```
> install.packages("renv")
> renv::restore()
```

# Data files and variables

[describe each column (variable) in each of your data files]

1. [data file name]
[list of variable anmes and descriptions]
2. [data file name]
[list of variable anmes and descriptions]
This data repository consist of xxxx data files, yyy code scripts, and this README document, with the following data and code filenames and variables


# Workflow

The shell scripts were executed on a 2021 iMac with Apple M1 chip 16GB memory and macOS version 14.7. 

The scripts included in this repository can reproduce the figures from raw reads, trait data, and map data.


## Visualization

The folder `plotting_scripts/` includes Rscripts each generates a figure or table in the manuscript.
To run all Rscripts, execute the following shell script

```
for file in plotting_scripts/*.R; do
    Rscript -e "renv::activate('.'); source('$file')"
done
```

Some figures use a pre-made cartoon (e.g., Figs 1, 2, 6, and S5), which is stored in `plots/cartoons/` as png/pdf format.


