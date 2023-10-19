#' This scripts remove the bottom 5% worst fastq reads via filtlong v0.2.1


library(tidyverse)
source(here::here("analysis/00-metadata.R"))

list_samples <- c()
i = 1
list.files(paste0(folder_data, "raw/Chang_Q5C_results/"))

paste0(folder_data, "raw/Chang_Q5C_results/Chang_Q5C_1/", "reads/raw_reads.fastq.gz")

# Create the folder structure
dir.create(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filter_reads/"), recursive = T)

# Invoke CLT
system(paste0(
    "cd; \n",
    "conda env list;",
    "filtlong --keep_percent 95 ", folder_data, "raw/Chang_Q5C_results/Chang_Q5C_1/reads/raw_reads.fastq.gz | gzip > ", folder_data, "temp/plasmidsaurus/Chang_Q5C_results/Chang_Q5C_1/01-filter_reads/output.fastq.gz"
))

