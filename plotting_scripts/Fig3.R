#' This script plots the core gene tree

renv::load()
library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
load(file = paste0(folder_data, "genomics_analysis/variants/snps.rdata"))


# core gene VA population ----
list_scaled_branches <- c(17, 18, 27, 13, 11, 12, 14)
pop = "VA"
p1 <- tr %>%
    drop.tip(isolates$genome_id[isolates$population != pop]) %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    geom_tiplab(align = T, hjust = -0.1, size = 3) +
    geom_tippoint(aes(color = site_group)) +
    #geom_label(aes(x = branch, label = node)) +
    #geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    #geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(20, -5, 0) *1e-4, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 8, width = 0.001) +
    scale_color_manual(values = site_group_colors) +
    scale_linetype_manual(values = c(1,5)) +
    scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = c(0.2, 0.9),
        legend.background = element_rect(fill = grey(0.9), color = NA),
        legend.key = element_blank(),
        legend.key.spacing.y = unit(0, "mm")
    ) +
    guides(linetype = "none", color = guide_legend(title = NULL) )+
    labs()

# core gene PA population ----
list_scaled_branches <- c(18, 31)
pop = "PA"
p2 <- tr %>%
    drop.tip(isolates$genome_id[isolates$population != pop]) %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    left_join(rename(isolates_contigs, label = genome_id)) %>%
    mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
    mutate(scaled_branch = ifelse(node %in% list_scaled_branches, T, F)) %>%
    mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
    as.treedata() %>%
    ggtree(aes(linetype = scaled_branch)) +
    geom_nodepoint(aes(label = highlight_boot), shape = 18, color = 1, size = 3, alpha = 0.3) +
    geom_tiplab(align = T, hjust = -0.1, size = 3) +
    geom_tippoint(aes(color = site_group)) +
    #geom_nodelab(aes(label = node)) +
    #geom_label(aes(x = branch, label = node)) +
    #geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("     ","     ","     "), label.r = unit(0.3, "lines"), label.size = 0, fill = c("#96BBBB", "#96BBBB", "#96BBBB"), alpha = 0.7) +
    #geom_label2(aes(subset=(node %in% c(31,37,39))), label = c("Ensifer spp.", "Ensifer meliloti", "Ensifer medicae"), label.size = 0, fill = NA, nudge_x = c(20, -5, 0) *1e-4, nudge_y = c(1, 1, -1), hjust = 1, fontface = "italic") +
    geom_treescale(x = 0, y = 8, width = 0.001) +
    scale_color_manual(values = site_group_colors) +
    scale_linetype_manual(values = c(1,5)) +
    scale_x_continuous(expand = c(0,0.001)) +
    theme_tree() +
    theme(
        legend.position = c(0.2, 0.9),
        legend.background = element_rect(fill = grey(0.9), color = NA),
        legend.key = element_blank(),
        legend.key.spacing.y = unit(0, "mm")
    ) +
    guides(linetype = "none", color = guide_legend(title = NULL)) +
    labs()

# SNPs ----
plot_snps <- function (isolates_i, pcoa_i, eigs_i) {
    isolates_i %>%
        bind_cols(tibble(mds1 = pcoa_i$points[,1], mds2 = pcoa_i$points[,2])) %>%
        ggplot() +
        geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
        geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
        geom_point(aes(x = mds1, y = mds2, color = site_group), size = 3, stroke = 1, shape = 21) +
        scale_color_manual(values = site_group_colors) +
        theme_classic() +
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            legend.position = "top",
            legend.title = element_blank()
        ) +
        guides(color = "none", shape = "none") +
        labs(x = paste0("PCoA axis 1(", eigs_i[1], "%)"), y = paste0("PCoA axis 1(", eigs_i[2], "%)"))

}


# Plotlist ----
plist <- list(
    p1,
    plot_snps(isolates1, pcoa1, eigs1),
    p2,
    plot_snps(isolates2, pcoa2, eigs2)
)


p <- plot_grid(plotlist = plist, nrow = 2, labels = LETTERS[c(1,3,2,4)],
               align = "h", axis = "lr",
               scale = 0.9, rel_widths = c(1.2, 1)) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig3.png"), p, width = 8, height = 6)
