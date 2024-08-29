#' This script plot the Fst per gene

renv::load()
library(tidyverse)
library(ggsci)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
sccg_fst <- read_csv(paste0(folder_data, "genomics_analysis/variants/sccg_fst.csv"))


p <- sccg_fst %>%
    filter(metric == "Gst_est") %>%
    select(-metric) %>%
    pivot_longer(-gene, names_to = "gradient") %>%
    mutate(gradient = ifelse(gradient == "elev", "elevation", "urbanization")) %>%
    ggplot() +
    geom_histogram(aes(x = value, fill = gradient), alpha = .5, position = "identity", binwidth = 0.01) +
    #scale_x_continuous(limits = c(-.05,1.05)) +
    scale_fill_aaas() +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(.8,.8),
        legend.background = element_rect(color = "grey10", fill = "white")
    ) +
    guides() +
    labs(x = "Nei's G_st", title = "988 core genes")
ggsave(paste0(folder_data, "genomics_analysis/variants/01-Fst_sccg.png"), p, width = 4, height = 4)

