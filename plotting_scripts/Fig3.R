#' GCV

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
library(lme4)
library(car)
source(here::here("metadata.R"))

load(paste0(folder_genomics, "pangenome/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
ani <- read_csv(paste0(folder_genomics, "taxonomy/ani.csv"))
isolates <- isolates %>% left_join(ani)
tt <- read_gpas()


# Panel A. core gene ----
tr <- tbtr$tr[[1]]
edges_to_scale <- which(tr$edge[,2] %in% 44)
tr$edge.length[edges_to_scale] <- tr$edge.length[edges_to_scale]*0.05

p1_1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    mutate(highlight = ifelse(node == 44, T, F)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse", aes(color = organism_name)) +
    geom_tippoint(aes(fill = organism_name), shape = 21, size = 1) +
    geom_nodepoint(aes(color = highlight, size = highlight)) +
    scale_color_manual(values = c(species_colors, `TRUE` = "red")) +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 0)) +
    scale_fill_manual(values = species_colors) +
    theme(
        plot.background = element_rect(color = "black", fill = "white", linewidth = 1)
    ) +
    guides(fill = "none", color = "none", size = "none")

p1 <- p1_1

# Panel B. gcv  ----
p2 <- tbtr$tr[[2]] %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = organism_name), hjust = 0, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = organism_name), shape = -1, size = -1) +
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 4, y = 25) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,3,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()


p2_1 <- isolates %>%
    filter(genome_id %in% get_taxa_name(p2)) %>%
    select(genome_id, region, organism_name) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
    ggplot() +
    geom_tile(aes(x = region, y = genome_id, fill = organism_name), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors, breaks = c("Sinorhizobium meliloti", "Sinorhizobium medicae")) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(5, "mm"),
        legend.text = element_text(face = "italic", size = 10),
        legend.box.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        plot.background = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2), nrow = 1)) +
    labs()

# Panel B. genes ----
tb <- tt$gpatl %>%
    filter(value == 1) %>%
    left_join(select(isolates, genome_id, organism_name))

genes_order <- tb %>%
    distinct(gene, genome_id) %>%
    group_by(gene) %>%
    count() %>%
    arrange(desc(n)) %>%
    pull(gene)

p3 <- tb %>%
    mutate(
        genome_id = factor(genome_id, rev(get_taxa_name(p2))),
        gene = factor(gene, genes_order),
        organism_name = factor(organism_name, c("Sinorhizobium meliloti", "Sinorhizobium medicae", "control"))
    ) %>%
    distinct(gene, organism_name, genome_id) %>%
    ggplot() +
    geom_tile(aes(x = gene, y = genome_id)) +
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradient(low = "grey80", high = "grey20") +
    facet_grid2(organism_name ~., scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "right",
        strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.background.x = element_rect(color = NA, fill = "gray90"),
        plot.background = element_blank()
    ) +
    guides() +
    labs(x = "gene family")

## Stat: do s meliloti and s medicate differ in their genes?
# mod <- tb %>%
#     filter(organism_name %in% c("S. medicae", "S. meliloti")) %>%
#     select(ge, genome_id, organism_name) %>%
#     group_by(ge, genome_id, organism_name) %>%
#     count() %>%
#     lmer(n ~ ge + organism_name + (1|organism_name:genome_id), data = .)
#Anova(mod, type = 3)
#emmeans(mod, ~organism_name)

# Panel C. functional genes -----
tb <- tt$gd %>%
    filter(str_detect(gene, "dna|grp|gro|rpo|clp|rec|uvr")) %>%
    mutate(ge = str_sub(gene, 1, 5) %>% str_remove("_\\d$|_$")) %>%
    mutate(g = str_sub(ge, 1, 3)) %>%
    select(g, ge, gene, genome_id) %>%
    left_join(select(isolates, genome_id, organism_name))

p4 <- tb %>%
    mutate(
        genome_id = factor(genome_id, rev(get_taxa_name(p2))),
        organism_name = factor(organism_name, c("Sinorhizobium meliloti", "Sinorhizobium medicae", "control"))
    ) %>%
    group_by(g, organism_name, genome_id, ge) %>%
    count() %>%
    mutate(n = factor(n)) %>%
    ggplot() +
    geom_tile(aes(x = ge, y = genome_id, fill = n)) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = setNames(grey(0.1*(7:1)), 1:7), name = "copy number") +
#    scale_fill_gradient(low = "grey80", high = "grey20", breaks = 1:10, name = "copy number") +
    facet_grid(organism_name ~ g, scales = "free", space = "free") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 0, size = 8, face = "italic"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        legend.position = "bottom",
        #legend.direction = "horizontal",
        strip.placement = "outside",
        strip.clip = "off",
        strip.text.x = element_text(face = "italic"),
        strip.text.y = element_blank(),
        strip.background.x = element_rect(color = NA, fill = "gray90")
    ) +
    guides(fill = guide_legend(nrow = 1, label.position = "bottom")) +
    labs()

# ----
p <- plot_grid(
    p2, p2_1 + guides(fill = "none"),
    p3, p4, nrow = 1,
    align = "h", axis = "tb", rel_widths = c(1,.2,1,2),
    labels = c("", "", "", "C")
) +
    draw_plot(p1, x = .01, y = .7, width = .13, height = .28) +
    draw_plot(get_legend(p2_1), x = -.35, y = -.42) +
    draw_text("A", x = .015, y = .96, size = 15, hjust = 0, fontface = "bold") +
    draw_text("B", x = .15, y = .96, size = 15, hjust = 0, fontface = "bold") +
    draw_text("core", x = .05, y = .95, size = 10, hjust = 0) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 6)

# PERMANOVA. Difference in heat relevant gene composition ?----
tbm <- tb %>%
    filter(organism_name %in% c("Sinorhizobium meliloti", "Sinorhizobium medicae")) %>%
    group_by(organism_name, genome_id, ge) %>%
    count() %>%
    pivot_wider(names_from = ge, values_from = n, values_fill = 0)
dm <- vegan::vegdist(tbm[,-c(1,2)], method = "bray")
vegan::adonis2(dm ~ organism_name, data = tbm, permutations = 999)
