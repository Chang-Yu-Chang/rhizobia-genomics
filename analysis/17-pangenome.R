#' This script analyzes the pangenome data

renv::load()
library(tidyverse)
library(janitor)
library(phangorn) # For NJ tree
source(here::here("analysis/00-metadata.R"))

pa <- read_delim(paste0(folder_genomics, "pangenome/panaroo/gene_presence_absence.Rtab"), delim = "\t")
m <- 1/tcrossprod(as.matrix(pa[,-1]))
tree <- phangorn::NJ(m)

#pa <- read_rtab(paste0(folder_genomics, "pangenome/panaroo/gene_presence_absence.Rtab"))
#g <- read_graph(paste0(folder_genomics, "pangenome/panaroo/final_graph.gml"), format = "gml")
#tree <- as.phylo(g)
#fit <- panstripe(pa, tree)
