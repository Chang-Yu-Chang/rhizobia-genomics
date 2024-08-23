
#' This script plots the Robinson-Foulds distance

renv::load()
library(tidyverse)
library(cowplot)
library(ape)
library(tidytree)
library(ggtree)
library(scales)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))

# Panel A trees ----
list_scaled_branches <- c(37:39, 14:16, 50)
scaling_factor = 20
pt1 <- tr_seq_core %>%
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
pt2 <- tr_gpa_genomes %>%
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
        #legend.background = element_rect(fill = grey(0.9), color = NA),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.spacing.y = unit(0, "mm"),
        legend.position = "top"
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs()

# Panel B. Normalized Robinson-Foulds tree distance ----
dist.topo(tr_seq_core, drop.tip(tr_gpa_genomes, "g20"), method = "PH85") %>% suppressWarnings()

# Random trees
trs <- rmtopology(1000, 30, tip.label = tr_seq_core$tip.label)
trs <- c(trs, tr_seq_core, drop.tip(tr_gpa_genomes, "g20"))
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
    coord_fixed(ratio = 1/8) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_log10(breaks = breaks_log(n = 5), labels = label_log()) +
    theme_classic() +
    theme(
        plot.margin = unit(c(0,0,0,0), "mm")
    ) +
    guides() +
    labs(x =  "Normalized Robinson-Foulds Distance", y = "Count")

p_cong <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], rel_widths = c(2,1), scale = c(0.95, 0.8)) +
    theme(plot.background = element_rect(color = NA, fill = "white"))


# Cophenetic permutation test ----
source(here::here("forposter/cophenetic.R"))
tbp <- tb1 %>%
    mutate(id = factor(id, 1:12)) %>%
    unnest(distances_pm) %>%
    compute_percentiles()
plot_cophenetic <- function (tbp, tb_obs) {
    tbp %>%
        ggplot() +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = alpha("gold", .2)) +
        geom_segment(aes(x = id, xend = id, y = p05, yend = p95), linewidth = 1,  arrow = arrow(length = unit(3, "mm"), angle = 90, ends = "both")) +
        geom_point(aes(x = id, y = p50, color = "95% CIs"), shape = 3, stroke = 1, size = 2) +
        geom_point(data = tb_obs, aes(x = id, y = distances_obs, color = "observation"), shape = 21, stroke = 2, size = 2) +
        scale_x_discrete(breaks = 1:12, labels = tb1$replicon_type) +
        scale_y_continuous(breaks = seq(0,1.5,0.5), limits = c(0, 1.5), minor_breaks = seq(0,1.5,0.1)) +
        scale_color_manual(values = c("observation" = "maroon", "95% CIs" = "black"), name = NULL) +
        facet_grid(.~ feature, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.8, 0.2),
            legend.background = element_rect(color = "black", fill = "white"),
            legend.box.background = element_rect(color = NA, fill = NA),
            legend.key = element_rect(color = NA, fill = NA),
            legend.spacing.y = unit(10,"mm"),
            legend.text = element_text(size = 8),
            legend.key.size = unit(5, "mm"),
            strip.background = element_blank(),
            strip.text = element_text(angle = 0, hjust = 0, vjust = 0, size = 8),
            strip.clip = "off",
            panel.spacing.x = unit(1, "mm"),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(size = 8, angle = 20, hjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 15),
            plot.margin = unit(c(10,10,5,5), "mm"),
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = "", y = expression(bar(d)["within"]/bar(d)["between"]))
}

p1 <- tbp %>%
    filter(population == "VA") %>%
    plot_cophenetic(filter(tb_obs, population == "VA")) +
    guides(color = "none") +
    labs(title = "Elevation gradient")
p2 <- tbp %>%
    filter(population == "PA") %>%
    plot_cophenetic(filter(tb_obs, population == "PA")) +
    labs(title = "Urbanization gradient")

p_coph <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb", label_y = .9, labels = c("C", "D"))

# Combine ----
p <- plot_grid(p_cong, p_coph, ncol = 1) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig5.png"), p, width = 10, height = 8)

