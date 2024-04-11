#' This script cleans the gene presence/absence and compute distance

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

gpa <- read_delim(paste0(folder_data, "genomics/pangenome/isolates/gene_presence_absence.Rtab"))
gpa <- clean_names(gpa)
write_csv(gpa, paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))

# Transpose the gene presence-absence table
gpat <- gpa %>%
    pivot_longer(cols = -gene, names_to = "genome_id") %>%
    pivot_wider(names_from = gene, values_from = value, values_fill = 0)

write_csv(gpat, paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
