#' This script calls the snp variation on chromosome core genomes
library(tidyverse)
library(cowplot)
library(janitor)
library(FactoMineR) # for MCA
library(ggsci)
source(here::here("analysis/00-metadata.R"))

# 0. read data ----
egc <- read_delim(paste0(folder_data, "/temp/anvio/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
