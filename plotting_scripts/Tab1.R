#' This script generates the blast results

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

ft <- iso %>%
    # Column cleanup
    select(...1, population, rrna_species, rrna_pident, contig_species, contig_pident, exp_id, genome_id, exp_lup, growth_curve) %>%
    mutate(population = ifelse(population == "VA", "Virginia", "Pennsylvania") %>% factor(c("Virginia", "Pennsylvania"))) %>%
    arrange(population, contig_species) %>%
    mutate(...1 = 1:n()) %>%
    flextable() %>%
    # Column labels
    set_header_labels(
        `...1` = "",
        exp_id = "Internal ID",
        population = "Region",
        genome_id = "Genome ID",
        contig_species = "Species", rrna_species = "Species",
        rrna_pident = "Identity (%)", contig_pident = "Identity (%)",
        exp_lup = "Plant inoculation",
        growth_curve = "Growth"
    ) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(j = c("rrna_species", "contig_species"), pr_t = fp_text_default(italic = T)) %>%
    add_header_row(values = c("", "BLAST rRNA gene", "BLAST genome", "Strain", "Experiment"), colwidths = c(2,2,2,2,2)) %>%
    # Align and spacing
    merge_v(j = c("population")) %>%
    valign(valign = "top") %>%
    align(align = "center", part = "all") %>%
    autofit() %>%
    # Font color
    style(j = "contig_species", pr_t = fp_text_default(bold = T, italic = T)) %>%
    color(color = species_colors["S. adhaerens"], j = "contig_species", i = ~contig_species == "S. adhaerens") %>%
    color(color = species_colors["S. canadensis"], j = "contig_species", i = ~contig_species == "S. canadensis") %>%
    color(color = "#32AEA0", j = "contig_species", i = ~contig_species == "S. medicae") %>%
    color(color = "#8C5488", j = "contig_species", i = ~contig_species == "S. meliloti") %>%
    # Lines and background
    hline(i = 15, j = 2:10) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = seq(1, nrow_part(.), 2), j = 3:10) %>%
    fix_border_issues()

#save_as_html(ft, path = here::here("plots/Tab1.html"))
save_as_image(ft, path = here::here("plots/Tab1.png"), res = 300)

