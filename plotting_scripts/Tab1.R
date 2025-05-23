#' This script generates the blast results

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))

ft <- iso %>%
    select(-exp_sativa) %>%
    flextable() %>%
    # Labels
    set_header_labels(
        `...1` = "",
        exp_id = "Strain",
        #gradient = "Gradient",
        population = "Population",
        genome_id = "Genome ID",
        contig_species = "Species", rrna_species = "Species",
        rrna_pident = "Identity (%)", contig_pident = "Identity (%)",
        exp_lup = "Plant inculation",
        growth_curve = "Growth"
    ) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(j = c("rrna_species", "contig_species"), pr_t = fp_text_default(italic = T)) %>%
    add_header_row(values = c("", "BLAST rRNA gene", "BLAST genome", "Experiment"), colwidths = c(4,2,2,2)) %>%
    # Align and spacing
    merge_v(j = c("population")) %>%
    valign(valign = "top") %>%
    align(align = "center", part = "all") %>%
    autofit() %>%
    # Lines and background
    hline(i = c(6,15,27), j = 2:10) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = seq(1, nrow_part(.), 2), j = 3:10) %>%
    fix_border_issues()


#save_as_html(ft, path = here::here("plots/Tab1.html"))
save_as_image(ft, path = here::here("plots/Tab1.png"), res = 300)

