#' This script generates the blast results

renv::load()
library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

iso <- isolates %>%
    left_join(isolates_tax) %>%
    # Remove these two lines after correcting the code
    mutate(contig_species = str_replace(contig_species, "E. ", "S. ")) %>%
    mutate(rrna_species = str_replace(rrna_species, "E. ", "S. ")) %>%
    mutate(population = ifelse(population == "VA", "elevation", "urbanization")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(population, site_group, genome_id) %>%
    select(population, site_group, genome_id, starts_with("rrna"), starts_with("contig")) %>%
    rename(site = site_group, genome = genome_id) %>%
    mutate(contig_length = round(contig_length/1000, 2),
           contig_pident = round(contig_pident, 1),
           rrna_pident = round(rrna_pident, 1)) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything()) %>%
    select(` `, population, site, genome, rrna_species, rrna_pident, contig_species, contig_pident)

# Flextable
ft <- flextable(iso) %>%
    valign(j = c("population", "site"), valign = "top") %>%
    set_header_labels(
        `...1` = "",
        population = "Gradient",
        site = "Site",
        genome = "Genome",
        contig_species = "Species", rrna_species = "Species",
        rrna_pident = "Identity (%)", contig_pident = "Identity (%)"
    ) %>%
    style(j = c("rrna_species", "contig_species"), pr_t = fp_text_default(italic = T)) %>%
    add_header_row(values = c("", "BLAST rRNA gene", "BLAST genome"), colwidths = c(4,2,2)) %>%
    bg(bg = "grey90", i = seq(1,38,2), j = 4:8) %>%
    align(align = "center", part = "all") %>%
    valign(valign = "center") %>%
    autofit() %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    merge_v(j = c("population", "site")) %>%
    hline(i = c(6,28), j = 3:8) %>%
    hline(i = c(16), j = 2:8) %>%
    fix_border_issues(part = "all") %>%
    autofit()

save_as_html(ft, path = here::here("plots/TabS7.html"))
