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
    mutate(population = ifelse(population == "VA", "mountain", "city")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(genome_id) %>%
    select(population, site_group, genome_id, starts_with("rrna"), starts_with("contig")) %>%
    rename(site = site_group, genome = genome_id) %>%
    mutate(contig_length = round(contig_length/1000, 2),
           rrna_pident = round(rrna_pident, 1),
           contig_pident = round(contig_pident, 1)) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

# Flextable
ft <- flextable(iso) %>%
    merge_v(j = c("population", "site")) %>%
    valign(j = c("population", "site"), valign = "top") %>%
    #separate_header() %>%
    set_header_labels(
        population = "Population",
        site = "Site",
        genome = "Genome",
        rrna_species = "Species", contig_species = "Species",
        rrna_sseqid = "Accession", contig_sseqid = "Accession",
        rrna_pident = "Idenity (%)", contig_pident = "Idenity (%)",
        rrna_length = "Length (bp)", contig_length = "Length (kbp)"
    ) %>%
    style(j = c("rrna_species", "contig_species"), pr_t = fp_text_default(italic = T)) %>%
    add_header_row(values = c("", "rRNA gene", "Genome"), colwidths = c(4, 4, 4)) %>%
    align(i = 1, j = NULL, align = "center", part = "header") %>%
    autofit() %>%
    hline_bottom() %>%
    fix_border_issues(part = "all")

save_as_image(ft, path = here::here("plots/TabS2.png"))

