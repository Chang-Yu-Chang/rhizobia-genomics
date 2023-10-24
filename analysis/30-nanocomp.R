#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data
nanostats <- read_delim(paste0(folder_data, "temp/plasmidsaurus/NanoStats.txt"), delim = ":", col_names = T)

# Clean the data
colnames(nanostats) <- c("summ", "gs")
nanostats <- nanostats %>% filter(!is.na(gs))

gs1 <- nanostats$gs[1:8] %>%
    str_split(" ") %>%
    lapply(function(x){tibble(genome_id = paste0("g", 1:19), value = x[x!=""])})

gs2 <- nanostats$gs[9:23] %>%
    str_split("\t") %>%
    lapply(function(x){tibble(genome_id = paste0("g", 1:19), value = x[x!=""])})

nanostats <- nanostats %>%
    mutate(summ_group = c(rep("s1", 8), rep("s2", 15))) %>%
    mutate(value = c(gs1, gs2)) %>%
    select(summ_group, summ, value) %>%
    unnest(cols = c(value))

# 1. plot the summary statisticcs
p <- nanostats %>%
    filter(summ_group == "s1") %>%
    mutate(value = as.numeric(str_replace_all(value, ",", ""))) %>%
    ggplot() +
    geom_histogram(aes(x = value, group = genome_id), color = "black", fill = NA) +
    facet_wrap(~summ, scales = "free_x", ncol = 2) +
    scale_x_continuous(labels = scales::label_comma(), n.breaks = 4) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/30-01-nanocomp_summary_statistics.png"), p, width = 6, height = 10)























