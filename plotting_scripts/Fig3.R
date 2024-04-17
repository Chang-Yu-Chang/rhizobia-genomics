#' This script plots

renv::load()
library(tidyverse)
library(cowplot)
library(ape)
library(tidytree)
library(ggtree)
library(proxy) # For computing jaccard distance
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
#library(vegan) # for computing jaccard
source(here::here("metadata.R"))

gpat <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpat.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

# Panel A-B. PAP heatmaps
ii <- isolates %>%
    left_join(isolates_contigs) %>%
    filter(!genome_id %in% c("g20", "g28")) %>%
    select(genome_id, site_group, population, species)

plot_heatmap <- function (pop = "VA", sgc = site_group_colors[1:2]) {
    gpatl %>%
        left_join(ii) %>%
        mutate(genome_id = factor(genome_id, rev(ii$genome_id))) %>%
        filter(population == pop) %>%
        ggplot() +
        geom_rect(data = tibble(site_group = names(sgc)), aes(fill = site_group), xmin = 0, xmax = 30000, ymin = 0, ymax = 100, alpha = 0.2) +
        geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = site_group_colors) +
        facet_grid(site_group~., scales = "free_y", space = "free_y") +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 8, color = "black"),
            strip.background = element_rect(color = NA, fill = NA),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
        ) +
        guides(fill = "none") +
        labs(x = "gene cluster", y = "genome", title = "Gene presence/absence")
}

p1 <- plot_heatmap()
p2 <- plot_heatmap("PA", site_group_colors[3:4])

# Panel C-D. matched trees
plot_matchedtree <- function (pop = "VA") {
    pop <- ifelse(pop == "VA", "PA", "VA")

    pt1 <- tr %>%
        drop.tip(isolates$genome_id[isolates$population == pop]) %>%
        as_tibble() %>%
        left_join(rename(ii, label = genome_id)) %>%
        as.treedata() %>%
        ggtree()

    # Plot accessory tree
    pt2 <- tr_acce %>%
        drop.tip(isolates$genome_id[isolates$population == pop]) %>%
        as_tibble() %>%
        left_join(rename(ii, label = genome_id)) %>%
        as.treedata() %>%
        ggtree()

    d1 <- pt1$data
    d2 <- pt2$data
    d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
    dd <- bind_rows(d1, d2) %>% filter(isTip)

    pt1 + geom_tree(data = d2) +
        geom_line(data = dd, aes(x, y, group = label), color = "grey90", linetype = 1, linewidth = .2) +
        geom_tippoint(aes(color = site_group), size = 2) +
        geom_tippoint(data = d2, aes(color = site_group), size = 2) +
        scale_x_continuous(limits = c(-1, 3)) +
        annotate("text", x = c(-0.5, max(pt2$data$x*1.5)), y = rep(max(pt1$data$y), 2), label = c("core genes", "gene content"), hjust = c(0.5,-0.5)) +
        scale_color_manual(values = site_group_colors) +
        theme_tree() +
        theme(
            legend.background = element_rect(color = "black", linewidth = .5)
        ) +
        guides(color = "none") +
        labs()

}

p3 <- plot_matchedtree()
p4 <- plot_matchedtree("PA")

# Panel E-F. SNPs
load(file = paste0(folder_data, "genomics_analysis/variants/snps.rdata"))

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
p5 <- plot_snps(isolates1, pcoa1, eigs1)
p6 <- plot_snps(isolates2, pcoa2, eigs2)


p <- plot_grid(p1, p3, p5, p2, p4, p6, nrow = 2, labels = LETTERS[c(1,3,5,2,4,6)],
               align = "h", axis = "lr",
               scale = 0.9, rel_widths = c(1, 1.5, 1)) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 6)




