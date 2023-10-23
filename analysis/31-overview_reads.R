#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data
list_g <- paste0(rep("Chang_Q5C_results", 19), "/Chang_Q5C_", 1:19, "/")
list_g[11] <- "Chang_Q5C_results_repeated/Chang_Q5C_11/"
list_g[18] <- "Chang_Q5C_results_repeated/Chang_Q5C_18/"
list_reads <- rep(list(NA), 19)

for (i in 1:19) list_reads[[i]] <- read_table(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "01-filtlong/raw_reads.txt"), col_names = c("name", "phred", "length"))
raw_reads <- list_reads %>%
    bind_rows(.id = "genome_id") %>%
    mutate(genome_id = factor(genome_id, 1:20))

# subtract by 33
raw_reads$phred <- raw_reads$phred-33


# 1. read length vs. phred
p <- raw_reads %>%
    ggplot() +
    geom_point(aes(x = length, y = phred), size = 0.2, alpha = 0.2) +
    scale_x_log10() +
    facet_wrap(~genome_id, ncol = 5, nrow = 4) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/31-01-read_length_vs_phred.png"), plot = p, width = 16, height = 20)

# 2. Summary stat for each file number histogram
raw_reads %>%
    group_by(genome_id) %>%
    summarize(n = n(),
              median_phred = median(phred),
              median_length = median(length))

# 3. plot read length for all reads
p <- raw_reads %>%
    filter(genome_id==1) %>%
    ggplot() +
    geom_histogram(aes(x = length)) +
    #scale_x_log10() +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/31-03-read_length.png"), plot = p, width = 4, height = 4)

# 4. plot read length for reads below 20k
p <- raw_reads %>%
    filter(length < 20000) %>%
    filter(genome_id==1) %>%
    ggplot() +
    geom_histogram(aes(x = length)) +
    #scale_x_log10() +
    #scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/31-04-read_length_20k.png"), plot = p, width = 4, height = 4)



























