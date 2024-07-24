#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))

# Gene content matrix ----
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))
gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpacl.csv")) %>%
    separate(contig_id, into = c("genome_id", "temp1", "temp2"), remove = F) %>%
    select(-temp1, -temp2) %>%
    mutate(gene = factor(gene, gene_order$gene)) %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))


p_gpam <- gpacl %>%
    ggplot() +
    geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
    scale_y_discrete(expand = c(0,0)) +
    facet_grid(~replicon_type, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        strip.background = element_blank(),
        strip.text = element_text(angle = 30, hjust = 0, vjust = 0),
        strip.clip = "off",
        plot.background = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")

#ggsave(paste0(folder_data, "genomics_analysis/gene_content/01-gpa_heatmap.png"), p, width = 6, height = 3)
nrow(gene_order) # 26886 genes in the pangenome


# Core gene tree ----
plot_tree <- function (tr, color_breaks = c("suburban", "urban", "bootstrap>95%")) {
    tr %>%
        as_tibble() %>%
        left_join(rename(isolates_contigs, label = genome_id)) %>%
        left_join(rename(isolates, label = genome_id)) %>%
        mutate(branch.length = ifelse(node %in% list_scaled_branches, branch.length * 0.01, branch.length)) %>%
        mutate(scaled_branch = ifelse(node %in% list_scaled_branches, "scaled to 1%", "not scaled")) %>%
        mutate(highlight_boot = ifelse(label > 95, T, F)) %>%
        as.treedata() %>%
        ggtree(aes(linetype = scaled_branch), linewidth = 1) +
        geom_tippoint(aes(color = site_group), size = 3) +
        geom_nodepoint(aes(label = highlight_boot, color = "bootstrap>95%"), shape = 16, size = 5, alpha = 0.3) +
        scale_color_manual(values = c(site_group_colors, `bootstrap>95%`="grey40"), name = NULL, breaks = color_breaks) +
        scale_linetype_manual(values = c(1,3), name = NULL) +
        # #coord_cartesian(clip = "off") +
        coord_flip(clip = "off") +
        theme_tree() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.8, 0.20),
            legend.background = element_rect(color = NA, fill = "grey90"),
            legend.box.background = element_rect(color = NA, fill = NA),
            legend.key = element_rect(color = NA, fill = NA),
            legend.spacing.y = unit(-3,"mm"),
            legend.text = element_text(size = 10),
            plot.margin = unit(c(0,10,0,0), "mm"),
            plot.background = element_blank(),
            panel.background = element_blank()
        ) +
        guides(linetype = "none") +
        labs()
}

list_scaled_branches <- c(15,16,17,18,27,13,12,11,14)
p1 <- tr_seq_core %>%
    drop.tip(isolates$genome_id[isolates$population == "PA"]) %>%
    plot_tree(c("high elevation", "low elevation"))

list_scaled_branches <- c(18, 31)
p2 <- tr_seq_core %>%
    drop.tip(isolates$genome_id[isolates$population == "VA"]) %>%
    plot_tree(c("suburban", "urban"))

p_core <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb", labels = c("A", "B"))

# Gene content tree ----
list_scaled_branches <- c(1,2,15,17,18,19)
p1 <- tr_gpa_genomes %>%
    drop.tip(isolates$genome_id[isolates$population == "PA"]) %>%
    plot_tree(c("high elevation", "low elevation")) +
    guides(color = "none")
    # geom_nodelab(aes(label = node)) +
    # geom_tiplab(aes(label = node))

list_scaled_branches <- c(18)
p2 <- tr_gpa_genomes %>%
    drop.tip(isolates$genome_id[isolates$population == "VA"]) %>%
    plot_tree(c("suburban", "urban")) +
    guides(color = "none")

p_gpa <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb", labels = c("C", "D"))



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
            legend.spacing.y = unit(-3,"mm"),
            legend.text = element_text(size = 8),
            legend.key.size = unit(1, "mm"),
            strip.background = element_blank(),
            strip.text = element_text(angle = 0, hjust = 0, vjust = 0, size = 8),
            strip.clip = "off",
            panel.spacing.x = unit(1, "mm"),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text(size = 8, angle = 20, hjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 15),
            plot.margin = unit(c(0,10,5,5), "mm"),
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = "", y = expression(bar(d)["within"]/bar(d)["between"]))
}


p1 <- tbp %>%
    filter(population == "VA") %>%
    plot_cophenetic(filter(tb_obs, population == "VA")) +
    guides(color = "none")
p2 <- tbp %>%
    filter(population == "PA") %>%
    plot_cophenetic(filter(tb_obs, population == "PA"))

p_coph <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb", labels = c("E", "F"), label_y = .9)


# Combine figures ----

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig4.png"), scale = 1) +
    draw_plot(p_gpam, width = .35, height = .35, x = 0, y = .35) +
    draw_plot(p_core, width = .5, height = .26, x = .45, y = .66) +
    draw_plot(p_gpa, width = .5, height = .26, x = .45, y = .35) +
    draw_plot(p_coph, width = .5, height = .36, x = .45, y = -0.05) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig4.png"), p, width = 15, height = 8)

