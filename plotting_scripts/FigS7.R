#' This script plots the kmer networks

renv::load()
library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(ggforce)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
dist_genetics <- read_csv(paste0(folder_data, "genomics_analysis/dist_genetics.csv"))
load(file = paste0(folder_data, "phylogenomics_analysis/networks/networks.rdata"))

#
plot_hist <- function (pop = "VA", d) {
    dist_genetics %>%
        left_join(select(isolates, genome_id1 = genome_id, population1 = population)) %>%
        left_join(select(isolates, genome_id2 = genome_id, population2 = population)) %>%
        filter(population1 == pop, population2 == pop) %>%
        filter(genome_id1 != genome_id2) %>%
        mutate(across(starts_with("genome_id"), ordered)) %>%
        filter(genome_id1 < genome_id2) %>%
        ggplot() +
        geom_histogram(aes(x = {{d}}), fill = "white", color = "black") +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        guides() +
        labs()
}
plot_net <- function (pop = "VA", d, d_thr) {
    g %>%
        activate(nodes) %>%
        left_join(isolates) %>%
        filter(population == pop) %>%
        activate(edges) %>%
        mutate(bin = ifelse({{d}} < d_thr, paste0("d_kmer<", d_thr), "nope")) %>%
        ggraph(layout = "linear", circular = T) +
        #geom_mark_hull(aes(x, y, group = species, fill = species), concavity = 4, expand = unit(3, "mm"), alpha = 0.25) +
        geom_edge_arc(aes(colour = bin), alpha = 0.1) +
        geom_node_point(aes(fill = site_group), shape = 21, color = "black", size = 5) +
        scale_edge_color_manual(values = c("black", "NA")) +
        scale_fill_manual(values = site_group_colors) +
        theme_void() +
        theme(
            plot.background = element_rect(color = NA, fill = "white"),
            plot.margin = unit(c(1,1,1,1), "mm"),
            legend.position = "right"
        ) +
        guides(edge_colour = "none") +
        labs()
}

#
plist <- list(
    plot_hist(pop = "VA", d = d_kmer) + labs(x = "Jaccard Index", title = "k-mers"),
    plot_net(pop = "VA", d = d_kmer, d_thr = 0.8),
    plot_hist(pop = "PA", d = d_kmer) + labs(x = "Jaccard Index", title = "k-mers"),
    plot_net(pop = "PA", d = d_kmer, d_thr = 0.8),
    plot_hist(pop = "VA", d = d_ani) + labs(x = "ANI", title = "ANI"),
    plot_net(pop = "VA", d = d_ani, d_thr = 0.1),
    plot_hist(pop = "PA", d = d_ani) + labs(x = "ANI", title = "ANI"),
    plot_net(pop = "PA", d = d_ani, d_thr = 0.1)
)

p <- plot_grid(plotlist = plist, ncol = 2, align = "v", axis = "lr", labels = LETTERS[1:8], scale = 0.9) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS7.png"), p, width = 8, height = 8)





if (FALSE) {
    # Panel B. Plot kmer heatmap ----
    thr_kmer <- 0.85
    sort_pairs <- function (edges, upper = T) {
        tt <- edges %>%
            mutate(pair = 1:n()) %>%
            pivot_longer(cols = c(from, to)) %>%
            group_by(pair) %>%
            mutate(value = factor(value, nodes$genome_id)) %>%
            arrange(value)

        if (upper) {
            tt <- mutate(tt, name = c("from", "to"))
        } else tt <- mutate(tt, name = c("to", "from"))
        return(pivot_wider(tt, names_from = name, values_from = value))

    }
    p2 <- edges %>%
        mutate(to = factor(to, nodes$genome_id)) %>%
        mutate(from = factor(from, rev(nodes$genome_id))) %>%
        ggplot() +
        geom_tile(aes(x = to, y = from, fill = d_kmer)) +
        #scale_fill_viridis(alpha = 0.8, direction = -1, begin = 0, end = 1, breaks = seq(0,1,0.25)) +
        scale_fill_gradient2(low = "gold", mid = "palegreen4", high = "steelblue", midpoint = 0.8) +
        scale_x_discrete(position = "top", drop = F) +
        scale_y_discrete(drop = F) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            axis.title = element_blank(),
            axis.text = element_blank()
        ) +
        labs(title = "Jaccard distance")
}
