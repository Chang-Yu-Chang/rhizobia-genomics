#' This script makes table1

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
quast <- read_csv(paste0(folder_data, "genomics/qcs/quast.csv"))
busco <- read_csv(paste0(folder_data, "genomics/qcs/busco.csv"))
checkm <- read_csv(paste0(folder_data, "genomics/qcs/checkm.csv"))
ani <- read_csv(paste0(folder_data, "genomics/taxonomy/ani.csv"))

isolates <- select(isolates, genome_id, exp_id, population, site, growth_curve)
quast <- select(quast, genome_id, gc_percent, total_length_10000_bp, n50, l50)
busco <- pivot_wider(busco, id_cols = genome_id, names_from = taxlevel, values_from = completeness, names_prefix = "busco_")
checkm <- rename(checkm, genome_id = bin_id, checkm_completeness = completeness, checkm_contamination = contamination)
ani <- select(ani, genome_id, organism_name, ani)

# Main table
tb <- isolates %>%
    left_join(checkm) %>%
    left_join(quast) %>%
    left_join(busco) %>%
    left_join(ani) %>%
    mutate(
        region = recode(population, "VA" = "Virginia", "PA" = "Pennsylvania"),
        region = factor(region, levels = c("Virginia", "Pennsylvania")),
        growth_curve = ifelse(growth_curve == 1, "+", ""),
        gc_percent = round(gc_percent, 1),
        across(starts_with("busco_"), ~ round(.x * 100, 1)),
        ani = round(ani, 1)
    ) %>%
    arrange(region, organism_name, genome_id) %>%
    mutate(ID = row_number()) %>%
    select(
        ID, region, genome_id, growth_curve,
        checkm_completeness, checkm_contamination,
        gc_percent, total_length_10000_bp, n50, l50,
        busco_class, busco_order, busco_family, busco_genus,
        organism_name, ani
    )

# Flextable
ft <- tb %>%
    flextable() %>%
    # Header labels
    set_header_labels(
        ID = "",
        region = "Region",
        genome_id = "Strain",
        growth_curve = "Growth assay",
        checkm_completeness = "Completeness (%)",
        checkm_contamination = "Contamination (%)",
        gc_percent = "GC (%)",
        total_length_10000_bp = "Length (≥10kb)",
        n50 = "N50",
        l50 = "L50",
        busco_class = "Class (432)",
        busco_order = "Order (639)",
        busco_family = "Family (1041)",
        busco_genus = "Genus (3013)",
        organism_name = "Species",
        ani = "ANI (%)"
    ) %>%
    # One clean header grouping
    add_header_row(
        values = c(
            "",
            "Isolate Information",      # population, genome_id, exp_id, growth_curve
            "CheckM",     # completeness + contamination
            "Quast",      # gc, length, n50, l50
            "BUSCO",      # busco_* columns
            "Taxonomy"                  # organism + ani
        ),
        colwidths = c(1, 3, 2, 4, 4, 2)   # total = 17
    ) %>%
    # Header styling
    style(part = "header", pr_t = fp_text_default(bold = TRUE)) %>%
    style(j = "organism_name", pr_t = fp_text_default(italic = TRUE)) %>%
    border_inner_h(part = "header", border = fp_border_default(width = 0)) %>%
    align(align = "center", part = "header") %>%
    align(align = "center", part = "all") %>%
    vline(j = c(4, 6, 10, 14), border = fp_border_default(color = "white", width = 3), part = "header") %>%
    # Backgrounds and borders
    bg(bg = "white", part = "body") %>%
    bg(bg = "grey95", i = 1, j = 1:16, part = "header") %>%
    bg(bg = "white", i = 2, j = 1:16, part = "header") %>%
    bg(bg = "grey95", i = seq(1, nrow_part(.), 2), j = 1:16)

save_as_html(ft, path = here::here("plots/Tab1.html"))
save_as_image(ft, path = here::here("plots/Tab1.png"), res = 300)
