#' This script plots the figures

library(tidyverse)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))


# Panel A. core gene ----
nodes_to_scale <- c(38, 40, 1, 2, 41, 42, 54)
tr <- tbtr$tr[[1]]
tr <- root(tr, outgroup = "g2")
edges_to_scale <- which(tr$edge[,2] %in% nodes_to_scale)
tr$edge.length[edges_to_scale] <- tr$edge.length[edges_to_scale]*0.01

p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    mutate(highlight = ifelse(node %in% nodes_to_scale, T, F)) %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    # geom_nodepoint(aes(fill = highlight), alpha = .5, size = 3, shape = 21, color = "white") +
    # scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    scale_x_continuous(limits = c(0, 0.008)) +
    geom_treescale(x = .004, y = 30) +
    facet_grid2(~` `) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    #theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"),
        legend.position = "inside",
        legend.position.inside = c(.2, .8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        plot.margin = unit(c(0,5,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()

# Panel B. gene content ----
p2 <- tbtr$tr[[2]] %>%
    as_tibble() %>%
    left_join(rename(iso, label = genome_id)) %>%
    mutate(` ` = "") %>%
    as.treedata() %>%
    ggtree(layout = "ellipse") +
    geom_tiplab(aes(label = label, color = contig_species), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
    geom_tippoint(aes(color = contig_species), shape = -1, size = -1) +
    scale_x_continuous(limits = c(0, 130)) +
    geom_treescale(x = 65, y = 30) +
    facet_grid2(~` `) +
    coord_cartesian(clip = "off") +

    theme_tree() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8),
        plot.margin = unit(c(0,5,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs()

# Panel C. gene content ----
tt <- read_gpas()
p3 <- tt$gpacl %>%
    left_join(isolates) %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
    mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
    drop_na(replicon_type) %>%
    mutate(rep_gene = paste0(replicon_type, gene)) %>%
    mutate(rep_gene = factor(rep_gene, paste0(tt$gene_order$replicon_type, tt$gene_order$gene))) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
    ggplot() +
    geom_tile(aes(x = rep_gene, y = genome_id, fill = contig_species)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors) +
    facet_grid2(~replicon_type, scales = "free", space = "free_x", switch = "y", strip = strip_vanilla(clip = "off")) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.spacing.x = unit(1, "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        #panel.background = element_rect(color = NA, fill = "black"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")


# ----
p <- plot_grid(
    p1, p2, p3, nrow = 1,
    scale = .9, rel_widths = c(1,1,1.5),
    labels = LETTERS[1:3]
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 5 )
