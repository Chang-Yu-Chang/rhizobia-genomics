#' This script plots the figures

library(tidyverse)
library(janitor)
library(cowplot)
library(ggh4x)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))


plot_tree <- function (tr_core) {
    tr_core %>%
        as_tibble() %>%
        left_join(rename(isolates, label = genome_id)) %>%
        mutate(` ` = "") %>%
        as.treedata() %>%
        ggtree() +
        #geom_nodepoint(color = "grey70", shape = 16, alpha = .5, size = 5) +
        geom_tiplab(aes(label = label, color = population, fill = population), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
        geom_tippoint(aes(color = population, fill = population), shape = -1, size = -1) +
        scale_color_manual(values = population_colors) +
        scale_fill_manual(values = population_colors) +
        scale_y_continuous(limits = c(1, length(tr_core$tip.label)), expand = c(0,.5)) +
        coord_cartesian(clip = "off") +
        facet_grid2(~` `) +
        theme_tree() +
        theme(
            # legend.position = "inside",
            legend.title = element_blank(),
            legend.background = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 8),
            plot.margin = unit(c(0,5,0,0), "mm")
        ) +
        guides(fill = guide_legend(override.aes = list(label = "", color = NA, size = 2)), color = guide_legend(override.aes = list(size = 2, shape = 21))) +
        labs()
    #labs(title = paste0("Elevation S. medicae\n", nrow(tt$list_sccg), " single-copy core genes"))
}
plot_heatmap <- function (tt, p_tree) {

    tt$gpacl %>%
        left_join(isolates) %>%
        mutate(genome_id = factor(genome_id, isolates$genome_id)) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        drop_na(replicon_type) %>%
        mutate(rep_gene = paste0(replicon_type, gene)) %>%
        mutate(rep_gene = factor(rep_gene, paste0(tt$gene_order$replicon_type, tt$gene_order$gene))) %>%
        mutate(genome_id = factor(genome_id, rev(get_taxa_name(p_tree)))) %>%
        ggplot() +
        geom_tile(aes(x = rep_gene, y = genome_id, fill = population)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_manual(values = population_colors) +
        facet_grid2(~replicon_type, scales = "free", space = "free_x", switch = "y") +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 10),
            strip.placement = "outside",
            panel.spacing.x = unit(1, "mm"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = unit(c(0,0,0,0), "mm")
        ) +
        guides(fill = "none") +
        labs(x = "gene cluster", y = "genome")
    #ggtitle(paste0("Total: ", n_all, ", Core: ", n_core, ", Accessory: ", n_accessory))
}
plot_ngen <- function (tt) {
    #tt<-read_gpas("elev_med")
    tb_ngenomes <- tt$gpacl %>%
        group_by(replicon_type, gene) %>%
        count(name = "n_genomes") %>%
        ungroup() %>%
        arrange(gene) %>%
        ungroup()

    tb_ngenomes %>%
        group_by(replicon_type, gene) %>%
        mutate(n_genomes_forplot = map(n_genomes, ~1:.x)) %>%
        unnest(n_genomes_forplot) %>%
        ungroup() %>%
        mutate(rep_gene = paste0(replicon_type, gene)) %>%
        mutate(rep_gene = factor(rep_gene, paste0(tt$gene_order$replicon_type, tt$gene_order$gene))) %>%
        mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce"))) %>%
        drop_na(replicon_type) %>%
        mutate(`shared genomes` = "shared genomes") %>%
        ggplot() +
        geom_tile(aes(x = rep_gene, y = n_genomes_forplot), fill = "black") +
        scale_y_continuous(breaks = c(5, 10, 15, 20), expand = c(0,0)) +
        facet_grid2(`shared genomes`~replicon_type, scales = "free", space = "free_x", strip = strip_vanilla(clip = "off"), margins = F) +
        coord_cartesian(clip = "off") +
        theme_classic() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            #strip.text.y = element_text(angle = 45),
            strip.text.y = element_blank(),
            strip.clip = "off",
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.length.x = unit(0, "mm"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
            panel.spacing.x = unit(1, "mm"),
            plot.background = element_blank(),
            plot.margin = unit(c(0,0,0,0), "mm"),
        ) +
        guides() +
        labs(x = "gene cluster", y = "# of genomes")
}


p_tree1 <- plot_tree(tr_elev_med_core) + geom_treescale(width = 1e-3, x = 0, y = 9)
p_tree2 <- plot_tree(tr_urbn_mel_core) + geom_treescale(width = 1e-3, x = 0, y = 16)
leg1 <- get_legend(p_tree1 + theme(legend.position = "right", legend.direction = "horizontal"))
leg2 <- get_legend(p_tree2 + theme(legend.position = "right", legend.direction = "horizontal"))
p_gpa1 <- plot_tree(tr_elev_med_gpa) + geom_treescale(width = 5, x = 0, y = 9)
p_gpa2 <- plot_tree(tr_urbn_mel_gpa) + geom_treescale(width = 5, x = 0, y = 16)
p_heat1 <- plot_heatmap(read_gpas("elev_med"), p_gpa1)
p_heat2 <- plot_heatmap(read_gpas("urbn_mel"), p_gpa2)
p_ngen1 <- plot_ngen(read_gpas("elev_med"))
p_ngen2 <- plot_ngen(read_gpas("urbn_mel"))

p <- plot_grid(
    p_tree1 + theme(legend.position = "none"),
    p_gpa1 + theme(legend.position = "none"),
    p_heat1,
    leg1, NULL, p_ngen1,
    NULL, NULL, NULL,
    p_tree2 + theme(legend.position = "none"),
    p_gpa2 + theme(legend.position = "none"),
    p_heat2,
    leg2, NULL, p_ngen2,
    ncol = 3,  align = "hv", axis = "blr", scale = .95,
    labels = c("A", "C", "E", rep("", 6), "B", "D", "F", rep("", 3)),
    rel_widths = c(1, 1, 2), rel_heights = c(10, 2, .5, 17, 2)
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig5.png"), p, width = 10, height = 6)
