#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
library(ggh4x)
library(ape)
library(tidytree)
library(ggtree)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

plot_tree <- function (tr) {
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
        labs()
    #labs(title = paste0("Elevation S. medicae\n", nrow(tt$list_sccg), " single-copy core genes"))
}

tbtr <- tbtr %>%
    mutate(p_tr = map2(tr, tr_type, ~plot_tree(.x) + ggtitle(.y)))

tbtr$p_tr[[1]] <- tbtr$p_tr[[1]] + theme(legend.position = "bottom")
tbtr$p_tr[[5]] <- tbtr$p_tr[[5]] + theme(legend.position = "bottom")

p <- plot_grid(
    plotlist = tbtr$p_tr, nrow = 2, scale = .95, align = "h", axis = "tb",
    labels = c("A", "", "", "", "B"), rel_heights = c(10, 17), rel_widths = c(1,1,1,1.3)
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/FigS7.png"), p, width = 10, height = 6)

