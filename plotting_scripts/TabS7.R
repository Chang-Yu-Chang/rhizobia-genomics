#' This script generates the blast results

renv::load()
library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
#iso <- read_csv(paste0(folder_data, "output/iso.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))

iso <- isolates %>%
    left_join(isolates_tax) %>%
    mutate(population = ifelse(population == "VA", "elevation", "urbanization")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    #arrange(genome_id) %>%
    arrange(population, site_group, genome_id) %>%
    select(population, site_group, genome_id, starts_with("sm"), starts_with("rrna"), starts_with("contig")) %>%
    rename(site = site_group, genome = genome_id) %>%
    mutate(contig_length = round(contig_length/1000, 2),
           sm_query_containment_ani = round(sm_query_containment_ani*100, 1),
           contig_pident = round(contig_pident, 1),
           rrna_pident = round(rrna_pident, 1)) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything()) %>%
    select(` `, population, site, genome, sm_species, sm_query_containment_ani, rrna_species, rrna_pident, contig_species, contig_pident)

# Flextable
ft <- flextable(iso) %>%
    valign(j = c("population", "site"), valign = "top") %>%
    set_header_labels(
        `...1` = "",
        population = "Gradient",
        site = "Site",
        genome = "Genome",
        sm_species = "Species", contig_species = "Species", rrna_species = "Species",
        sm_query_containment_ani = "Query ANI (%)", rrna_pident = "Idenity (%)", contig_pident = "Idenity (%)"
    ) %>%
    style(j = c("sm_species", "rrna_species", "contig_species"), pr_t = fp_text_default(italic = T)) %>%
    add_header_row(values = c("", "Sourmash", "BLAST rRNA gene", "BLAST genome"), colwidths = c(4,2,2,2)) %>%
    bg(bg = "grey90", i = seq(1,38,2), j = 4:10) %>%
    align(align = "center", part = "all") %>%
    valign(valign = "center") %>%
    autofit() %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    delete_columns(j = c("sm_species", "sm_query_containment_ani")) %>%
    merge_v(j = c("population", "site")) %>%
    hline(i = c(6,16,28), j = 3:10) %>%
    fix_border_issues(part = "all")

save_as_html(ft, path = here::here("plots/TabS7.html"))
