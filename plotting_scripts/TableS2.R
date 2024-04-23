#' This script generates Table S2

renv::load()
library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

iso <- isolates %>%
    left_join(isolates_tax) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    select(population, site_group, genome_id, starts_with("rrna"), starts_with("contig")) %>%
    mutate(hit_rrna = paste0(rrna_species, " (", rrna_sseqid, ") ", round(rrna_pident, 1), "%, ", rrna_length)) %>%
    mutate(hit_contig = paste0(contig_species, " (", contig_sseqid, ") ", round(contig_pident, 1), "%, ", contig_length)) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything(), -starts_with("rrna"), -starts_with("contig"))

# Flextable
ft <- flextable(iso) %>%
    autofit() %>%
    merge_v(j = c("population", "site_group")) %>%
    valign(j = c("population", "site_group"), valign = "top")
    #style(j = 2, pr_t = fp_text_default(italic = T))

save_as_image(ft, path = here::here("plots/TabS2.png"))
