#' This script generates the blast results

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
isolates_all <- read_csv(paste0(folder_data, "mapping/isolates_all.csv"))

iso <- isolates %>%
    filter(!genome_id %in% c("g28", "g38")) %>% # duplicated of g20 and g40
    left_join(isolates_tax) %>%
    # Remove these two lines after correcting the code
    mutate(contig_species = str_replace(contig_species, "E. ", "S. ")) %>%
    mutate(rrna_species = str_replace(rrna_species, "E. ", "S. ")) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    arrange(gradient, population, genome_id) %>%
    select(gradient, population, exp_id, genome_id, starts_with("rrna"), starts_with("contig")) %>%
    mutate(contig_length = round(contig_length/1000, 2),
           contig_pident = round(contig_pident, 1),
           rrna_pident = round(rrna_pident, 1)) %>%
    mutate(` ` = 1:n()) %>%
    left_join(select(isolates_all, exp_id, genome_id, exp_lup, exp_sativa, growth_curve)) %>%
    mutate(
        exp_lup = ifelse(exp_lup == 1, "+", ""),
        exp_sativa = ifelse(exp_sativa == 1, "+", ""),
        growth_curve = ifelse(growth_curve == 1, "+", ""),
    ) %>%
    select(` `, population, exp_id, genome_id, rrna_species, rrna_pident, contig_species, contig_pident, exp_sativa, exp_lup, growth_curve)

write_csv(iso, paste0(folder_data, "output/iso.csv"))


ft <- flextable(iso) %>%
    # Labels
    set_header_labels(
        `...1` = "",
        exp_id = "Strain",
        #gradient = "Gradient",
        population = "Population",
        genome_id = "Genome ID",
        contig_species = "Species", rrna_species = "Species",
        rrna_pident = "Identity (%)", contig_pident = "Identity (%)",
        exp_sativa = "Alternate host",
        exp_lup = "Source host",
        growth_curve = "Growth"
    ) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(j = c("rrna_species", "contig_species"), pr_t = fp_text_default(italic = T)) %>%
    add_header_row(values = c("", "BLAST rRNA gene", "BLAST genome", "Experiment"), colwidths = c(4,2,2,3)) %>%
    # Align and spacing
    merge_v(j = c("population")) %>%
    valign(valign = "top") %>%
    align(align = "center", part = "all") %>%
    autofit() %>%
    # Lines and background
    hline(i = c(6,15,27), j = 2:11) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = seq(1, nrow_part(.), 2), j = 3:11) %>%
    fix_border_issues()


save_as_html(ft, path = here::here("plots/TabS1.html"))
save_as_image(ft, path = here::here("plots/TabS1.png"), res = 300)

