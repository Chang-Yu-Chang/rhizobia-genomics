#' This script make the distance matrix for NJ trees

renv::load()
library(tidyverse)
library(janitor)
library(ape)
library(ggtree)
library(tidytree)
source(here::here("analysis/00-metadata.R"))


gd <- read_csv(paste0(folder_data, "genomics/pangenome/panaroo/gene_data.csv"))

gd %>%
    distinct(clustering_id)
