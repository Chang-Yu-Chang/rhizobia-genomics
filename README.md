# Rhizobia genomics

Scripts for genome assembly, annotation, pangenome analysis, and phylogenies For Chang et al. 2024 manuscript


# Setup 

This repository relies heavily on `conda`/`mamba` for managing bioinformatic tools and on `renv` for R packages. 

### conda/mamba

- Follow commands in `setup_conda.sh` to install conda/miniforge/mamba on your terminal.

- Once you have installed `mamba`, follow commands in `setup_envs.sh` to install bioconda packages. In general, each bioconda tool has one mamba virtual environment.

### renv

The Rscripts are run under the following R version

```
> sessionInfo()
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.4
```

For first time `renv` user, install renv and restore the packages recorded in `renv.lock`

```
> install.packages("renv")
> renv::restore()
```
