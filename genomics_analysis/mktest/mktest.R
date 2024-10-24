#' This script performs MK test on each gene

library(tidyverse)
library(iMKT)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))



set_name = "elev_med"
#set_name = "urbn_mel"
# tt <- read_gpas(set_name)
# ff <- read_fsts(set_name)
#
# ff$per_gene_fst %>% view

gene = "zwf"
#gene = "accA"
gene = "IMPDH_3"
daf <- read_delim(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/", gene, ".daf"))
div <- read_delim(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/tables/", gene, ".div"))
standardMKT(daf, div)

# myDafData
# myDivergenceData
# standardMKT(myDafData, myDivergenceData)

