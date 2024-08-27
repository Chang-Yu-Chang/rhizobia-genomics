#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
#contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))

# Gene content matrix ----
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))
gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpacl.csv")) %>%
    separate(contig_id, into = c("genome_id", "temp1", "temp2"), remove = F) %>%
    select(-temp1, -temp2)
background_df <- tibble(site_group = c("high elevation", "low elevation", "suburban", "urban"))

p_gpam <- gpacl %>%
    filter(genome_id != "g42") %>%
    drop_na(replicon_type) %>%
    left_join(isolates) %>%
    mutate(gene = factor(gene, gene_order$gene)) %>%
    mutate(genome_id = factor(genome_id, rev(isolates$genome_id))) %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = site_group), alpha = 0.2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_tile(aes(x = gene, y = genome_id), fill = "grey10") +
    scale_fill_manual(values = site_group_colors) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0), position = "right") +
    facet_grid(site_group ~ replicon_type, scales = "free", space = "free", switch = "y" ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 15),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        strip.background = element_blank(),
        panel.spacing.y = unit(c(0,3,0), "mm"),
        strip.text.x = element_text(angle = 15, hjust = 0, vjust = 0, size = 15),
        #strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0),
        strip.text.y.left = element_text(angle = 75, hjust = .8, vjust = .5, size = 15),
        strip.clip = "off",
        plot.background = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "gene cluster", y = "genome")

nrow(gene_order) # 26504 genes in the pangenome

## How many genes are shared in chrome, pA and pB
gpacl %>%
    drop_na(replicon_type) %>%
    group_by(replicon_type, gene) %>%
    count() %>%
    mutate(is_core = n == 36) %>%
    group_by(replicon_type, is_core) %>%
    count

gpacl %>%
    #drop_na(replicon_type) %>%
    group_by(gene) %>%
    count() %>%
    mutate(is_core = n == 36) %>%
    group_by(is_core) %>%
    count



# Core gene tree ----
plot_tree <- function (tr, color_breaks = c("suburban", "urban", "bootstrap>95%")) {
    tr %>%
        as_tibble() %>%
        #left_join(rename(isolates_contigs, label = genome_id)) %>%
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
            #legend.position = "inside",
            legend.position = "none",
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

p_core <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb")

# Gene content tree ----
list_scaled_branches <- c(1,2,15,17,18,19)
p1 <- tr_gpa_genomes %>%
    drop.tip(isolates$genome_id[isolates$population == "PA"]) %>%
    plot_tree(c("high elevation", "low elevation")) +
    guides(color = "none")

list_scaled_branches <- c(18)
p2 <- tr_gpa_genomes %>%
    drop.tip(isolates$genome_id[isolates$population == "VA"]) %>%
    plot_tree(c("suburban", "urban")) +
    guides(color = "none")

p_gpa <- plot_grid(p1, p2, nrow = 1, scale = .95, align = "h", axis = "tb")


#  ----
p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig4.png"), scale = 1) +
    # draw_plot(p_gpam, width = .65, height = .3, x = 0.3, y = .05) +
    draw_plot(p_gpam, x = 0.05, y = .02, width = .35, height = .54) +
    draw_plot(p_core, x = .48, y = .5, width = .5, height = .38) +
    draw_plot(p_gpa, x = .48, y = .05, width = .5, height = .38) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig4.png"), p, width = 15, height = 8)

