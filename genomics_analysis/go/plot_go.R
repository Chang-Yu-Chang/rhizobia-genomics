#' This script plots the GO enrichment analysis

library(tidyverse)
library(janitor)
library(ggsci)
library(cowplot)
source(here::here("metadata.R"))

read_goenrich <- function (set_name) {
    goenrich_bygene <- read_csv(paste0(folder_data, "genomics_analysis/go/", set_name, "/goenrich_bygene.csv")) %>%
        mutate(goidterm = paste0(go_id, " ", term))
    goenrich_bysnp <- read_csv(paste0(folder_data, "genomics_analysis/go/", set_name, "/goenrich_bysnp.csv")) %>%
        mutate(goidterm = paste0(go_id, " ", term))

    return(list(goenrich_bygene = goenrich_bygene, goenrich_bysnp = goenrich_bysnp))
}
plot_goenrich <- function (goenrich) {
    p1 <- goenrich %>%
        mutate(goidterm = factor(goidterm, rev(goenrich$goidterm))) %>%
        ggplot() +
        geom_col(aes(x = goidterm, y = annotated, fill = category)) +
        coord_flip() +
        scale_fill_aaas() +
        theme_bw() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(.9,.9),
            legend.background = element_rect(color = "black"),
            axis.text.y.left = element_text(hjust = 0)
        ) +
        guides() +
        labs(x = "", title = set_name)

    p2 <- goenrich %>%
        mutate(goidterm = factor(goidterm, rev(goenrich$goidterm))) %>%
        ggplot() +
        geom_col(aes(x = goidterm, y = -log(as.numeric(classic_fisher), 10), fill = category)) +
        geom_hline(yintercept = -log(0.05, 10)) +
        coord_flip() +
        scale_fill_aaas() +
        theme_bw() +
        theme(
            axis.text.y = element_blank(),
            legend.position = "none"
        ) +
        guides() +
        labs(x = "", y = "-log(p)")
    p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, .8))
    return(p)
}

set_name = "elev_med"
ee <- read_goenrich(set_name)
p <- plot_goenrich(ee$goenrich_bygene)
ggsave(paste0(folder_data, "genomics_analysis/go/", set_name,"-01-go_bygene.png"), p, width = 15, height = 8)
p <- plot_goenrich(ee$goenrich_bysnp)
ggsave(paste0(folder_data, "genomics_analysis/go/", set_name,"-02-go_bysno.png"), p, width = 15, height = 8)

set_name = "urbn_mel"
ee <- read_goenrich(set_name)
p <- plot_goenrich(ee$goenrich_bygene)
ggsave(paste0(folder_data, "genomics_analysis/go/", set_name,"-01-go_bygene.png"), p, width = 15, height = 8)
p <- plot_goenrich(ee$goenrich_bysnp)
ggsave(paste0(folder_data, "genomics_analysis/go/", set_name,"-02-go_bysno.png"), p, width = 15, height = 8)
