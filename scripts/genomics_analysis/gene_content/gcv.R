#'

library(tidyverse)
library(cowplot)
library(vegan)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
symbiosis_genes <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/symbiosis_genes.csv"))
tt <- read_gpas()

# Symbiosis gene content similarity ----
# Test if co-occurring genomes share similar symbiosis genome
gpa_sym <- tt$gd %>%
    filter(str_remove(gene, "_\\d+") %in% symbiosis_genes$gene) %>%
    select(genome_id, gene) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = genome_id, values_fill = 0) %>%
    arrange(gene)

s1 <- iso$genome_id[iso$population == "PA" & iso$contig_species == "S. meliloti"]
s2 <- iso$genome_id[iso$population == "VA" & iso$contig_species == "S. meliloti"]
s3 <- iso$genome_id[iso$population == "PA" & iso$contig_species == "S. medicae"]
s4 <- iso$genome_id[iso$population == "VA" & iso$contig_species == "S. medicae"]

# Within sp, between populations
dist_longer <- function (d) {
    d %>%
        as.matrix() %>%
        as.tibble() %>%
        mutate(row = names(.)) %>%
        pivot_longer(cols = -row, names_to = "col", values_to = "d") %>%
        mutate(row = ordered(row, iso$genome_id), col = ordered(col, iso$genome_id)) %>%
        filter(row != col) %>%
        filter(row > col)
}
tbd <- vegdist(t(gpa_sym[,-1]), method = "bray", diag = F) %>%
    dist_longer()

# within sp, within pop
d1 <- c(
    tbd$d[tbd$row %in% s1 & tbd$col %in% s1],
    tbd$d[tbd$row %in% s2 & tbd$col %in% s2],
    tbd$d[tbd$row %in% s3 & tbd$col %in% s3],
    tbd$d[tbd$row %in% s4 & tbd$col %in% s4]
)


# within sp, between pop
d2 <- c(
    tbd$d[tbd$row %in% s1 & tbd$col %in% s2],
    tbd$d[tbd$row %in% s3 & tbd$col %in% s4]
)

# between sp, within pop
d3 <- c(
    tbd$d[tbd$row %in% s1 & tbd$col %in% s3],
    tbd$d[tbd$row %in% s2 & tbd$col %in% s4]
)

# between sp, between pop
d4 <- c(
    tbd$d[tbd$row %in% s1 & tbd$col %in% s4],
    tbd$d[tbd$row %in% s2 & tbd$col %in% s3]
)

bind_rows(
    tibble(t = "1in in", d = d1),
    tibble(t = "2in bet", d = d2),
    tibble(t = "3bet in", d = d3),
    tibble(t = "4bet bet", d = d4)
) %>%
    ggplot() +
    geom_boxplot(aes(x = t, y = d)) +
    geom_jitter(aes(x = t, y = d))
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()

t.test(d3, d4)

#
p1 <- tbtr$tr[[3]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = 0, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 4, y = 5) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()

#
tb <- tt$gd %>%
    #filter(!str_detect(gene, "group")) %>%
    filter(str_remove(gene, "_\\d+") %in% symbiosis_genes$gene) %>%
    #filter(str_detect(gene, "exo|rkp|lpt|rgt|nod|nif|fix")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(iso, genome_id, contig_species))

p2 <- tb %>%
    mutate(
        genome_id = factor(genome_id, rev(get_taxa_name(p1))),
        contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
    ) %>%
    group_by(g, contig_species, genome_id, ge) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = ge, y = genome_id, fill = n)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "grey80", high = "grey20") +
    #facet_grid2(contig_species ~ g, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "right",
        strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.background.x = element_rect(color = NA, fill = "gray90")
    ) +
    guides() +
    labs()

d <- dist(t(gpa_sym[,-1]))
tr <- hclust(d, method = "ward.D2") %>% as.phylo()
plot(tr)

p2

p <- plot_grid(
    p1, p2,
    align = "h", axis = "tb"
)


# Plot heat relevant genes
tb <- tt$gd %>%
    filter(!str_detect(gene, "group")) %>%
    #filter(str_remove(gene, "_\\d+") %in% symbiosis_genes$gene) %>%
    filter(str_detect(gene, "dna|grp|gro|rpo|clp|rec|uvr")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(iso, genome_id, contig_species))

tb %>%
    mutate(
        genome_id = factor(genome_id, rev(get_taxa_name(p1))),
        contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens", "control"))
    ) %>%
    group_by(g, contig_species, genome_id, ge) %>%
    count() %>%
    ggplot() +
    geom_tile(aes(x = ge, y = genome_id, fill = n)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "grey80", high = "grey20") +
    #facet_grid2(contig_species ~ g, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "right",
        strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.background.x = element_rect(color = NA, fill = "gray90")
    ) +
    guides() +
    labs()
