#' This script plot the Fst per gene and per snp

library(tidyverse)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_gcv_fsts <- function (set_name) {
    per_acce_fst <- read_csv(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"/per_acce_fst.csv"))
    per_genome_fst <- read_csv(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"/per_genome_fst.csv"))
    return(list(per_acce_fst = per_acce_fst, per_genome_fst = per_genome_fst))
}
set_name <- "elev_med"
ff <- read_gcv_fsts(set_name)
tt <- read_gpas(set_name)
gene_replicon <- tt$gpacl %>%
    select(gene, replicon_type) %>%
    replace_na(list(replicon_type = "others")) %>%
    distinct()
gene_order <- unique(tt$gpacl$gene)


acce_fst <- ff$per_acce_fst %>%
    drop_na(Gst) %>%
    mutate(gene = factor(gene, gene_order)) %>%
    left_join(gene_replicon) %>% # get the replicon where the gene is
    group_by(replicon_type) %>%
    mutate(replicon_type = factor(replicon_type, c("chromosome", "pSymA", "pSymB", "pAcce", "others"))) %>%
    ungroup()

acce_n <- acce_fst %>%
    filter(replicon_type %in% c("chromosome", "pSymA", "pSymB", "pAcce")) %>%
    group_by(replicon_type) %>%
    count()


p <- acce_fst %>%
    #filter(Fst == "Gprime_st") %>%
    filter(replicon_type %in% c("chromosome", "pSymA", "pSymB", "pAcce")) %>%
    ggplot() +
    geom_point(aes(x = gene, y = Gprime_st), shape = 21, size = .5) +
    #geom_text(data = acce_n, aes(label = paste0(n ," genes")), x = Inf, y = -Inf, vjust = -1, hjust = 1.2) +
    facet_grid2(~replicon_type, space = "free_x", scales = "free_x", strip = strip_vanilla(clip = "off")) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
    ) +
    guides() +
    labs(x = "gene cluster", y = "G'st")
ggsave(paste0(folder_data, "genomics_analysis/gcv_fst/", set_name,"-01-acce_fst_gprime.png"), p, width = 8, height = 3)


