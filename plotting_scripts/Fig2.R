#' This script plots the figures

library(tidyverse)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
library(gggenomes)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
seqs <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/seqs.csv"))
gens <- read_csv(paste0(folder_data, "phylogenomics_analysis/comparative/gens.csv"))


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
    scale_color_manual(values = species_colors) +
    scale_x_continuous(limits = c(0, 0.008)) +
    geom_treescale(x = .004, y = 30) +
    facet_grid2(~` `) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"),
        legend.position = "inside",
        legend.position.inside = c(.2, .8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides(color = "none") +
    labs()

p1_1 <- isolates %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    select(genome_id, population, contig_species) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p1)))) %>%
    ggplot() +
    geom_tile(aes(x = population, y = genome_id, fill = contig_species), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
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
    scale_color_manual(values = species_colors) +
    geom_treescale(x = 65, y = 30) +
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

p2_1 <- isolates %>%
    left_join(select(iso, genome_id, contig_species)) %>%
    select(genome_id, population, contig_species) %>%
    mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
    ggplot() +
    geom_tile(aes(x = population, y = genome_id, fill = contig_species), color = "black", linewidth = .5) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_manual(values = species_colors) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.clip = "off",
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        panel.background = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,-1), "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(linewidth = .2))) +
    labs()

# Panel C. pangenome ----
p3 <- gggenomes(
    seqs = mutate(seqs, bin_id = factor(bin_id, get_taxa_name(p2))) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "others"))),
    genes = filter(gens, str_detect(name, "nif|nod|fix"))
) +
    geom_bin_label() +
    geom_seq(aes(color = replicon_type), linewidth = 1) +
    geom_gene(aes(fill = name)) +
    scale_color_manual(values = c(chromosome = "black", pSymA = "darkred", pSymB = "darkblue", others = "grey80"), name = "Replicon") +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    theme_tree() +
    theme(
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.position = "top",
        legend.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()

# ----
p <- plot_grid(
    p1, p1_1 + guides(fill = "none"),
    p2, p2_1 + guides(fill = "none"),
    p3, nrow = 1,
    scale = .95, rel_widths = c(1,.2,1,.2,1.5),
    align = "h", axis = "tb",
    labels = c("A", "", "B", "C")
) +
    draw_text("Single-copy core gene", x = .05, y = .9, size = 10, hjust = 0) +
    draw_text("Gene content variation", x = .35, y = .9, size = 10, hjust = 0) +
    draw_plot(get_legend(p1_1), x = -.45, y = .25) +
    theme(plot.background = element_rect(color = NA, fill = "white"))


ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 6)


#
nrow(tt$gpa) # 26544

core <- tt$gpatl %>%
    group_by(gene) %>%
    filter(value == 1) %>%
    count() %>%
    ungroup() %>%
    filter(n == max(n))
nrow(core) / nrow(tt$gpa) *100

tt$gpacl %>%
    filter(str_detect(gene, "nod")) %>%
    filter(genome_id %in% paste0("g", c(2,3,15)))




if (F) {
    # Panel C. gene content ----

    tt <- read_gpas()
    p3 <- tt$gpacl %>%
        left_join(isolates) %>%
        left_join(select(iso, genome_id, contig_species)) %>%
        #mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
        #mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        #drop_na(replicon_type) %>%
        #mutate(rep_gene = paste0(replicon_type, gene)) %>%
        #mutate(rep_gene = factor(rep_gene, paste0(tt$gene_order$replicon_type, tt$gene_order$gene))) %>%
        mutate(gene = factor(gene, unique(tt$gene_order$gene))) %>%
        mutate(genome_id = factor(genome_id, rev(get_taxa_name(p2)))) %>%
        ggplot() +
        geom_tile(aes(x = gene, y = genome_id, fill = contig_species)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = species_colors) +
        #facet_grid2(~replicon_type, scales = "free", space = "free_x", switch = "y", strip = strip_vanilla(clip = "off")) +
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
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            #axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0,0,0,-1), "mm")
        ) +
        guides(fill = "none") +
        labs(x = "Cluster of orthologous genes", y = "genome")


}
