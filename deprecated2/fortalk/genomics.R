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
        #geom_tiplab(aes(label = species), offset = 1e-4) +
        geom_tippoint(aes(color = site_group), size = 3) +
        geom_nodepoint(aes(label = highlight_boot, color = "bootstrap>95%"), shape = 16, size = 5, alpha = 0.3) +
        #geom_nodelab(aes(label = node)) +
        #geom_tiplab(aes(label = node)) +
        geom_treescale() +
        scale_color_manual(values = c(site_group_colors, `bootstrap>95%`="grey40"), name = NULL, breaks = color_breaks) +
        scale_linetype_manual(values = c(1,3), name = NULL) +
        coord_cartesian(clip = "off") +
        #coord_flip(clip = "off") +
        theme_tree() +
        theme(
            #legend.position = "inside",
            legend.position = "top",
            legend.position.inside = c(0.8, 0.20),
            legend.background = element_blank(),
            legend.box.background = element_rect(color = NA, fill = NA),
            legend.key = element_rect(color = NA, fill = NA),
            legend.spacing.y = unit(-3,"mm"),
            legend.text = element_text(size = 10),
            plot.margin = unit(c(0,10,0,0), "mm"),
            plot.background = element_rect(color = NA, fill = "white"),
            panel.background = element_blank()
        ) +
        guides(linetype = "none") +
        labs()
}

list_scaled_branches <- c(15,16,17,18,27,13,12,11,14)
p_core <- tr_seq_core %>%
    drop.tip(isolates$genome_id[isolates$population == "PA"]) %>%
    drop.tip("g20") %>%
    plot_tree(c("high elevation", "low elevation"))


# Gene content tree ----
#list_scaled_branches <- c(1,2,15,17,18,19)
list_scaled_branches <- c(1, 2, 3, 14, 17, 18, 19)
#list_scaled_branches <- c(0)
p_gpa <- tr_gpa_genomes %>%
    drop.tip(isolates$genome_id[isolates$population == "PA"]) %>%
    drop.tip("g20") %>%
    plot_tree(c("high elevation", "low elevation"))


ggsave(here::here("fortalk/tree1.png"), p_core, width = 5, height = 5)
ggsave(here::here("fortalk/tree2.png"), p_gpa, width = 5, height = 5)

