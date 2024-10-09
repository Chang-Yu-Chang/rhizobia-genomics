#' Concatenate the single copy core genes

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

# Elevation melicae
set_name <- "elev_med"
list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"))



gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv"))
gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv"))
gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))

gpatl %>%
    #group_by(gene, genome_id) %>%
    distinct(gene, genome_id) %>%
    group_by(gene) %>%
    count() %>%
    arrange(desc(n))

gpacl$annotation_id


ncol(gpa)
list_sccg <- gpa$gene[which(apply(gpa[,-1], 1, sum) == ncol(gpa)-1)]
length(list_sccg)
tb_sccg_fst <- tibble(gene = list_sccg, fst = NA)
folder_set <- paste0(folder_data, "genomics_analysis/fst/set1/")
