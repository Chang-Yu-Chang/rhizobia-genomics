# General

Scripts and minimal data for reproducing the anaylsis, tables, and figures presented in the manuscript

"Biogeographic and genomic signatures of thermal adaptation in facultative symbionts"

Chang, Chang-Yu, Terrence Topping-Brown, Jazmine L. Rud, Mccall B. Calvert, Gerardo Bencosme, and Corlett Wood. 

# Setup

This repository relies heavily on `conda`/`mamba` for managing bioinformatic tools and on `renv` for R packages.

### conda/mamba

- Follow commands in `setup_conda.sh` to install conda/miniforge/mamba on your terminal.

- Once you have installed `mamba`, follow commands in `setup_envs.sh` to install bioconda packages. In general, each bioconda tool has one mamba virtual environment.

### renv

The Rscripts are executed under the following R environment

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

For first time `renv` user, install renv and restore the packages recorded in `renv.lock`

```
> install.packages("renv")
> renv::restore()
```

# Data folders and files

There are data folders. The data included in this repository should allow a user to reproduce all the figures and tables in both main text and supplements.

Some folders are labeled with (NOT IN THIS REPO) because it exceeds the github repo file limit. These files will be available in the final version.

1. `raw/` 
    - `growth_curves/`: the raw growth assay data
    - `plants/`: sampling site coordinate and plant experiment data
    - `ensifer_ncbi.csv`: the list of Sinorhizobium/Ensifer reference genomes
    - (NOT IN THIS REPO) `plasmidsaurus/`: contains the raw ON long-read sequences from Plasmidsaurus
2. `mapping/` has a mapping file
    - `isolates.csv` is the list of isolates/strains and their sequencing batches and ids
3. `genomics/` contains the intermediate files of genome assembly, annotation, and pangenomes
4. `genomics_analysis/`
    - `contigs/contigs.csv`: contig data
    - `genomes/`
        - `genomes.csv`: genome size
        - `qcs.csv`: quality control metrics from quast and busco
    - `distances/distl.csv`: long-form ANI and kmer data
    - `taxonomy/`: intermediate files for taxonomic identification using blast
        - `isolates_tax.csv`: contig level and rRNA level blast result for strains
    - The folders below are dependent on panaroo pangenomic output, so each folder consists of two subforders
        - `elev_med/`: S medicae strains in the elevation gradient
        - `urbn_mel/`: S meliloti strains in the urbanization gradient
    - `gene_content/`: gene presence and absence. Basically cleaned panaroo tables
    - `fst/` and `gcv_fst/`: fst for SNPs and GCV (gene content variation)
    - `dxy/` and `gcv_dxy/`: dxy for SNPs and GCV
    - `go/` and `gcv_go`: GO enrichment analysis for SNPs and GCV
    - `trees/trees.rdata`: R phylo objects of whole-genome tree
    - `replicon_trees/trees.rdata`: R phylo objects of replicon-level tree
6. `phenotypes/` 
    - `sites/`: field sampling sites
    - `growth/`: growth traits


# Workflow

The shell scripts were executed on a 2021 iMac with Apple M1 chip 16GB memory.

## Reproducing figures and tables

The folder `plotting_scripts/` includes Rscripts each generates a figure or table in the manuscript. Some figures use a pre-made cartoon (e.g., Figs 1, 2, 6, and S5), which is stored in `plots/cartoons/` as png/pdf format.

To run all Rscripts, execute the following shell script:

```
for file in plotting_scripts/*.R; do
    Rscript -e "renv::activate('.'); source('$file')"
done
```

## Reproducing anaylsis from raw data

To reproduce analysis from raw data (raw reads and trait data), analysis, to the final figures and tables, execute the following shell script:

```
# Assembly
cd genomics/assembly/
zsh assess_reads.sh             # Filter and output the filtered read txt
zsh denovo_assembly.sh          # Assembly
zsh assess_assemblies.sh        # Quality control (quast and busco); the busco mamba env binary needs to be specified in zshrc
zsh manual_concat.sh            # Manually concatenate the two genomes g20 and g24
cd ../../
Rscript -e "renv::activate('.'); source('genomics/assembly/consolidate_qcs.R')" # Consolidates quality control metrics

# Annotation
cd genomics/annotation/
zsh annotate_genomes.sh
zsh consolidate_annotations.sh 

# Taxonomy
cd genomics/taxonomy
zsh download_ncbi_genomes.sh  # Download genomes from NCBI
zsh ani.sh # Compare the pairwise distance between our assembled genomes and the reference genomes
cd ../../
Rscript -e "renv::activate('.'); source('genomics/taxonomy/consolidate_ani.R')"

# Pangenome
cd genomics/pangenome/
zsh pangenome.sh 
cd ../../
Rscript -e "renv::activate('.'); source('genomics/pangenome/clean_gpa.R')"           # Clean up panaroo outputs, mostly gene presence/absence; output csv files are noted in the R script
cd genomics/pangenome/
zsh concatenate_alignment.sh # Concatenate the single-copy core genes
zsh compute_trees1.sh # Use iqtree to compute single-copy core-gene trees
cd ../../
Rscript -e "renv::activate('.'); 
source('genomics/pangenome/compute_trees2.R')   # computes trees for gene content variation
source('genomics/pangenome/curate_trees.R')     # combines the tree objects into a Rdata file `trees.rdata`
source('genomics/pangenome/rf_tree.R')          # computes the RF distance between trees



# Fst and dxy and go
Rscript -e "renv::activate('.'); 
source('genomics_analysis/fst/compute_fst.R');          # Compute Fst for SNPs in core genes
source('genomics_analysis/gcv_fst/compute_gcv_fst.R');  # Compute Fst for accessory gene content variation (presence/absence)
source('genomics_analysis/dxy/compute_dxy.R');          # Compute dxy for SNPs in core genes
source('genomics_analysis/gcv_dxy/compute_gcv_dxy.R');  # Compute dxy for accessory gene content variation (presence/absence)
source('genomics_analysis/go/go.R');                    # Perform GO analysis for core genes
source('genomics_analysis/gcv_go/gcv_go.R');            # Perform GO analysis for accessory genes
"

# Growth assay
Rscript -e "renv::activate('.'); source('phenotypes/growth/fit_gc.R')"  # Smooth the raw growth curve and computes the growth traits 

# Map and climate
Rscript -e "renv::activate('.'); source('phenotypes/sites/extract_climate.R')" # Use DAYMET https://daymet.ornl.gov/ database to extract the climate data for our sampling sites given the coordinates
```
