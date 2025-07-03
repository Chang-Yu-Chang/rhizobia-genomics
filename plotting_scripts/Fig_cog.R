#'

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
library(flextable)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
tt <- read_gpas()
tb_cog <- tibble(genome_id = factor(iso$genome_id, iso$genome_id)) %>%
    left_join(select(iso, genome_id, population, contig_species)) %>%
    mutate(cog_classify = map(genome_id, ~read_tsv(paste0(folder_data, "genomics/cog/", .x, "/cog_classify.tsv")))) %>%
    unnest(cog_classify) %>%
    clean_names()

#

tb_cog %>%
    group_by(genome_id, cog_letter) %>%
    count() %>%
    pivot_wider(names_from = cog_letter, values_from = n) %>%
    mutate(genome_id = factor(genome_id, get_taxa_name(p2))) %>%
    arrange(genome_id) %>%
    flextable() %>%
    mk_par(j = 2, value = as_paragraph(minibar(value = A, max = max(A))),part = "body") %>%
    mk_par(j = 3, value = as_paragraph(minibar(value = B, max = max(B))),part = "body") %>%
    mk_par(j = 4, value = as_paragraph(minibar(value = C, max = max(C))),part = "body") %>%
    mk_par(j = 5, value = as_paragraph(minibar(value = D, max = max(D))),part = "body") %>%
    mk_par(j = 6, value = as_paragraph(minibar(value = E, max = max(E))),part = "body") %>%
    mk_par(j = 7, value = as_paragraph(minibar(value = F, max = max(F))),part = "body") %>%
    mk_par(j = 8, value = as_paragraph(minibar(value = G, max = max(G))),part = "body") %>%
    mk_par(j = 9, value = as_paragraph(minibar(value = H, max = max(H))),part = "body") %>%
    mk_par(j = "X", value = as_paragraph(minibar(value = X, max = max(X))),part = "body")






# Panel C. cog ----
# p4 <- tb_cog %>%
#     mutate(
#         genome_id = factor(genome_id, rev(get_taxa_name(p2))),
#         contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
#     ) %>%
#     group_by(cog_letter, cog_id, contig_species, genome_id) %>%
#     #filter(cog_letter %in% LETTERS[20:26]) %>%
#     count() %>%
#     ggplot() +
#     geom_tile(aes(x = cog_id, y = genome_id), fill = "maroon", color = NA) +
#     coord_cartesian(clip = "off") +
#     facet_grid2(~cog_letter, scales = "free_x", space = "free_x") +
#     theme_bw() +
#     theme(
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.spacing.x = unit(0, "cm")
#     ) +
#     guides() +
#     labs()

