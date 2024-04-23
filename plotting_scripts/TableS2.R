#' This script generates Table S2

renv::load()
library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
b_16s <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/b_16s.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))

isolates %>%
    select(population, site_group, site, genome_id)

colnames(en) <- c("Accession", "Species", "Strain")
en <- en %>%
    mutate(Species = str_replace(Species, "E.", "Ensifer"))

# Flextable
ft <- flextable(en) %>%
    autofit() %>%
    style(j = 2, pr_t = fp_text_default(italic = T))

save_as_image(ft, path = here::here("plots/TabS1.png"))
