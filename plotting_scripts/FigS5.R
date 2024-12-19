#' This script plots replicon trees

library(tidyverse)
library(cowplot)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/replicon_trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
contigs <- read_csv(paste0(folder_data, "genomics_analysis/contigs/contigs.csv"))
n_sccg <- read_csv(paste0(folder_data, "phylogenomics_analysis/replicon_trees/n_sccg.csv"))
n_accessory <- read_csv(paste0(folder_data, "phylogenomics_analysis/replicon_trees/n_accessory.csv"))

n_sccg <- n_sccg %>% mutate(tr_type = "core") %>% rename(n_gene = n_sccg)
n_accessory <- n_accessory %>% mutate(tr_type = "gpa") %>% rename(n_gene = n_accessory) %>% select(-n_core)
n_genes = bind_rows(n_sccg, n_accessory)

plot_tree <- function (tr, n_gene) {
    tr %>%
        as_tibble() %>%
        left_join(rename(isolates, label = genome_id)) %>%
        mutate(` ` = "") %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(aes(label = label, color = population, fill = population), hjust = -.1, align = T, offset = 1e-3, linetype = 3, linesize = .1) +
        geom_tippoint(aes(color = population, fill = population), shape = -1, size = -1) +
        scale_color_manual(values = population_colors) +
        scale_fill_manual(values = population_colors) +
        scale_y_continuous(limits = c(1, length(tr$tip.label)), expand = c(0,.5)) +
        coord_cartesian(clip = "off") +
        theme_tree() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            legend.margin = margin(0,0,0,0,"mm"),
            legend.box.margin = margin(0,0,0,0,"mm"),
            legend.background = element_blank(),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = .2),
            plot.margin = margin(0,5,0,0, "mm")
        ) +
        guides(fill = guide_legend(override.aes = list(label = "", color = NA, size = 2)), color = guide_legend(override.aes = list(size = 2, shape = 21))) +
        labs(title = n_gene)
}

tbtr <- tbtr %>%
    left_join(n_genes) %>%
    mutate(p_tr = map2(tr, n_gene, plot_tree))

tbtr$p_tr[[1]] <- tbtr$p_tr[[1]] + theme(legend.position = "bottom")
tbtr$p_tr[[9]] <- tbtr$p_tr[[9]] + theme(legend.position = "bottom")


list_p <- rep(list(NULL), 16)
list_p[1:4] <- tbtr$p_tr[1:4]
list_p[5:7] <- tbtr$p_tr[9:11]
list_p[9:12] <- tbtr$p_tr[5:8]
list_p[13:16] <- tbtr$p_tr[12:15]

p_combined <- plot_grid(
    plotlist = list_p, ncol = 4, scale = .9, align = "h", axis = "tb",
    labels = c("A", "", "", "", "B", "", "", "", "C", "", "", "", "D", "", "", ""),
    rel_heights = c(10, 17, 10, 17), rel_widths = c(1,1,1,1), label_x = -.1
)

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/FigS5.png"), scale = 1) +
    draw_plot(p_combined, x = .1, y = 0, width = .9, height = .95) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/FigS5.png"), p, width = 10, height = 12)

