#'

library(tidyverse)
source(here::here("metadata.R"))
set.seed(42)

iso <- read_csv(paste0(folder_data, "output/iso.csv"))
tt <- read_gpas()

genome_ids <- c("g43", "g37", "g33", "g42", "g27", "g26", "g25", "g36", "g24", "g10", "g20", "g23", "g22", "g41", "g32", "g45", "g31", "g34", "g35", "g44", "g39", "g21", "g19", "g9", "g17", "g16", "g5", "g4", "g8", "g13", "g11", "g6", "g30", "g29", "g40", "g3", "g2", "g15")

p <- tt$gd %>%
    filter(!str_detect(gene, "group")) %>%
    filter(str_detect(gene, "nod|nif|fix")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    mutate(
        genome_id = factor(genome_id, rev(genome_ids)),
        contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
    ) %>%
    group_by(g, contig_species, genome_id, ge) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = ge, y = genome_id, fill = n)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "grey80", high = "grey20") +
    facet_grid(contig_species ~ g, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "top"
    ) +
    guides() +
    labs()

 ggsave(paste0(folder_data, "genomics_analysis/gene_content/sym.png"), p, width = 5, height = 5)



 tt$gd







tt$gpa %>%
    filter(!str_detect(gene, "group")) %>%
    filter(str_detect(gene, "gro")) %>%
    view


tt$gpatl %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    filter(!str_detect(gene, "group")) %>%
    filter(str_detect(gene, "gro|dnaK")) %>%
    filter(value == 1) %>%
    group_by(contig_species, genome_id) %>%
    count() %>%
    view



tt$gpatl %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    filter(!str_detect(gene, "group")) %>%
    filter(str_detect(gene, "dnaK")) %>%
    filter(value == 1) %>%
    group_by(contig_species, genome_id) %>%
    count() %>%
    view



