#' This script performs MK test

library(tidyverse)
library(iMKT)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_fsts <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}

set_name = "elev_med"
#set_name = "urbn_mel"
tt <- read_gpas(set_name)
ff <- read_fsts(set_name)

ff$per_gene_fst %>% view

myDafData
myDivergenceData
standardMKT(myDafData, myDivergenceData)

