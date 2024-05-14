#' This script plots the Robinson-Foulds distance

renv::load()
library(tidyverse)
library(cowplot)
library(ape)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))

# Panel A trees ----
list_scaled_branches <- c(37:39, 14:16, 50)
scaling_factor = 20
pt1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    left_join(rename(isolates, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    # scaling factor for matching tree x scales
    mutate(branch.length = branch.length * scaling_factor) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(-5, -5, 0) *1e-4*scaling_factor, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    scale_linetype_manual(values = c(1,5)) +
    scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(linetype = "none") +
    labs()

# Plot accessory tree
list_scaled_branches <- c(1,2,31,33,34,35,46)
pt2 <- tr_gpa %>%
    as_tibble() %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    left_join(rename(isolates, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    #geom_tiplab(align = T, size = 3) +
    geom_label2(aes(subset=(node %in% c(32,35,46))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(aes(subset=(node %in% c(32,35,46))), label = c("Ensifer spp.", "Ensifer medicae", "Ensifer meliloti"), label.size = 0, fill = NA, nudge_x = c(60, -10, -10)*1e-3, nudge_y = c(1, 1, 1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 28, width = 0.01) +
    scale_color_manual(values = species_colors) +
    scale_linetype_manual(values = c(1,5)) +
    theme_tree() +
    theme(
        legend.position = "top",
        plot.margin = unit(c(0,10,0,0), "mm")
    ) +
    guides(linetype = "none", color = "none") +
    labs()

d1 <- pt1$data
d2 <- pt2$data
## reverse x-axis and set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.05
dd <- bind_rows(d1, d2) %>% filter(isTip)

scale_width = 0.01
p1 <- pt1 +
    geom_tree(data = d2, aes(linetype = scaled_branch)) +
    ggnewscale::new_scale_fill() +
    geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
    # left tree
    geom_tippoint(aes(color = site_group), size = 2) +
    annotate("segment", x = -0.05, xend = -0.05 + scale_width, y = 25, yend = 25) +
    annotate("text", x = -0.05+scale_width/2, y = 25, label = format(scale_width / scaling_factor, scientific = F), vjust = -1, hjust = 0.5) +
    # right tree
    geom_tippoint(data = d2, aes(color = site_group), size = 2) +
    geom_nodepoint(data = d2, aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    geom_label2(data = d2, aes(subset=(node %in% c(32,35,46))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    geom_label2(data = d2, aes(subset=(node %in% c(32,35,46))), label = c("Ensifer spp.", "Ensifer medicae", "Ensifer meliloti"), label.size = 0, fill = NA, nudge_x = c(-10, 110, 100)*1e-3, nudge_y = c(1, 1, 1), hjust = 1, fontface = "italic") +
    annotate("segment", x = 0.5, xend = 0.5+scale_width, y = 25, yend = 25) +
    annotate("text", x = 0.5+scale_width/2, y = 25, label = scale_width, vjust = -1, hjust = 0.5) +
    scale_color_manual(values = site_group_colors) +
    scale_x_continuous(limits = c(-0.1, 0.55), expand = c(0,0)) +
    # Label
    annotate("text", x = c(-0.05, 0.5), y = c(30,30), label = c("core genes", "gene content")) +
    theme_tree() +
    theme(
        legend.background = element_rect(fill = grey(0.9), color = NA),
        legend.key = element_blank(),
        legend.key.spacing.y = unit(0, "mm"),
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs()

# Panel B. Normalized Robinson-Foulds tree distance ----
dist.topo(tr, drop.tip(tr_gpa, "g20"), method = "PH85") %>% suppressWarnings()

# Random trees
trs <- rmtopology(1000, 30, tip.label = tr$tip.label)
trs <- c(trs, tr, drop.tip(tr_gpa, "g20"))
md <- dist.topo(trs)
dim(as.matrix(md)) # 1002

tbb <- as_tibble(as.matrix(md)/max(md)[1]) %>%
    mutate(row = paste0("tree", 1:1002)) %>%
    pivot_longer(cols = -row, names_to = "col", values_drop_na = T)

obs <- tbb %>% filter(row == "tree1001", col == "tree1002") %>% pull(value)

p2 <- tbb %>%
    filter(row < col) %>%
    filter(row != "tree1001") %>%
    arrange(value) %>%
    ggplot() +
    geom_histogram(aes(x = value), breaks = seq(0, 1, 0.05), color = "black", fill = "white") +
    geom_vline(xintercept = obs, color = "maroon", linetype = 2) +
    scale_y_log10() +
    coord_fixed(ratio = 1/8) +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides() +
    labs(x =  "Normalized Robinson-Foulds Distance", y = "Count")

# Combine ----
p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], rel_widths = c(2,1), scale = c(0.95, 0.8)) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 4)

