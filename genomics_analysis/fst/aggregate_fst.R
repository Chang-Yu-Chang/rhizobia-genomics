#' This script computes Fst for a set of genes

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

set1_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/set1_fst.csv"))

set1_fst <- set1_fst %>%
    filter(metric == "Gst_est")
    # group_by(gene) %>%
    # nest()


set1_fst$pop_fst <- NA

for (i in 1:length(unique(set1_fst$gene))) {
    cat("\n", i)
    set1_fst$pop_fst[i] <- list(read_csv(paste0(folder_data, "genomics_analysis/fst/set1/", list_sccg[i], "-pop.csv"),  col_types = cols()))
}

set1_pop_fst <- set1_fst %>%
    unnest(pop_fst) %>%
    select(gene, singlecopy, pop1, pop2, Gst_est)

write_csv(set1_pop_fst, paste0(folder_data, "genomics_analysis/fst/set1_pop_fst.csv"))
