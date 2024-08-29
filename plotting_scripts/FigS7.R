#' This script plots the genomic distance matrices and mantel test

renv::load()
library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
dist_genomes <- read_csv(paste0(folder_data, "genomics_analysis/distances/dist_genomes.csv"))
tb_mantel <- read_csv(paste0(folder_data, "genomics_analysis/distances/tb_mantel.csv"))
dist_genomes <- dist_genomes %>%
    mutate(genome_id1 = ordered(genome_id1, genomes$genome_id), genome_id2 = ordered(genome_id2, genomes$genome_id)) %>%
    filter(genome_id1 >= genome_id2) %>%
    mutate(genome_id1 = ordered(genome_id1, rev(genomes$genome_id)), genome_id2 = ordered(genome_id2, genomes$genome_id))

# Plot distance matrices
plot_dist <- function(d_i) {
    dist_genomes %>%
        ggplot() +
        geom_tile(aes(x = genome_id2, genome_id1, fill = {{d_i}}), width = .8, height = .8) +
        geom_tile(aes(x = genome_id2, genome_id1), color = "black", fill = NA) +
        scale_x_discrete(position = "bottom", expand = c(0,0)) +
        scale_y_discrete(position = "left", expand = c(0,0)) +
        scale_fill_gradient(low = "snow", high = "darkred", name = "distance") +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "inside",
            legend.position.inside = c(.8, .8),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.margin = unit(c(0,0,0,10), "mm")
        ) +
        guides() +
        labs(x = "genome", y = "genome")
}

plist <- list(
    plot_dist(d_kmer) + ggtitle("k-mer (31-mer)"),
    plot_dist(d_sccg) + ggtitle("988 single-copy core genes"),
    plot_dist(d_ani) + ggtitle("1-ANI"),
    plot_dist(d_jaccard) + ggtitle("Jaccard distance on gene presence/absence")
)


# Plot the network for mantel test
p_mantel <- tb_mantel %>%
    mutate(corr = round(corr,2)) %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(node_label = str_replace(name, "d_", "")) %>%
    ggraph(layout = "linear", circular = T) +
    geom_edge_arc2(aes(edge_width = corr)) +
    scale_edge_width(range = c(.5,2), name = "Mantel's r") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
        plot.margin = unit(c(10,10,10,10), "mm"),
        panel.background = element_rect(fill = alpha("grey", 0.3), color = NA),
        legend.position = "inside",
        legend.position.inside = c(1,.5),
        legend.key.height = unit(0, "mm")
    ) +
    guides()


#
p_dist <- plot_grid(plotlist = plist, ncol = 2, align = "hv", axis = "tblr", labels = LETTERS[1:4], label_x = .2, label_y = .9, scale = 0.75) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

p <- ggdraw(p_dist) +
    draw_plot(p_mantel, x = .4, y = .4, width = .25, height = .25)

ggsave(here::here("plots/FigS7.png"), p, width = 10, height = 10)

