#' This script analysis the anvio generated

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 0. read data
egc <- read_delim(paste0(folder_data, "/temp/anvio/summary/ensifer_gene_clusters_summary.txt"), delim = "\t", show_col_types = F)
isolates_anvio <- read_delim(paste0(folder_data, "temp/42-isolates_anvio.txt"), show_col_types = F) %>%
    select(sample, r_category)
isolates <- read_csv(paste0(folder_data, "temp/42-isolates.csv"), show_col_types = F) %>%
    mutate(sample = str_replace(genome_id, "g", "Chang_Q5C_")) %>%
    select(sample, everything(), -r, -lag, -maxOD) %>%
    bind_rows(tibble(sample = c("em1021", "em1022", "wsm419"))) %>%
    left_join(isolates_anvio) %>%
    mutate(strain = ifelse(is.na(strain), "ncbi", strain),
           strain_site = ifelse(is.na(strain_site), "ncbi", strain_site),
           strain_site_group = ifelse(is.na(strain_site_group), "ncbi", strain_site_group))


egc %>%




#
egc %>%
    replace_na(list(bin_name = "rest")) %>%
    group_by(gene_cluster_id, bin_name) %>%
    summarize(count = n())

