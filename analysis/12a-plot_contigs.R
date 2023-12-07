#' This script summarize assembled genome information

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(janitor)
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F) %>%
    mutate(genome_id = factor(genome_id, genome_id))
contigs <- read_csv(paste0(folder_data, "temp/12-contigs.csv"), show_col_types = F) %>%
    left_join(isolates) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))
contigs_large <- read_csv(paste0(folder_data, "temp/12-contigs_large.csv"), show_col_types = F) %>%
    left_join(isolates) %>%
    mutate(genome_id = factor(genome_id, isolates$genome_id))


# 1. Number of contigs ----
p1 <- contigs %>%
    group_by(genome_id) %>%
    count(name = "n_contigs") %>%
    ggplot() +
    geom_histogram(aes(x = n_contigs), color = "black", fill = "white", binwidth = 1) +
    scale_x_continuous(breaks = seq(0, 100, 2), limits = c(1,40)) +
    scale_y_continuous(breaks = seq(0, 10, 1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "# of contigs", title = "all contigs")

p2 <- contigs_large %>%
    group_by(genome_id) %>%
    count(name = "n_contigs") %>%
    ggplot() +
    geom_histogram(aes(x = n_contigs), color = "black", fill = "white", binwidth = 1) +
    scale_x_continuous(breaks = seq(0, 100, 1), limits = c(1,10)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "# of contigs", title = "contigs > 500k")

p <- plot_grid(p1, p2, nrow = 2, align = "v", axis = "rl", labels = c("A", "B"))
ggsave(paste0(folder_data, "temp/12a-01-n_contigs.png"), plot = p, width = 4, height = 8)


# 2. plot the genome size by contigs ----

p1 <- contigs %>%
    group_by(genome_id) %>%
    mutate(contig_ordered = factor(1:n())) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = contig_length/10^6, group = contig_ordered, fill = contig_ordered), position = position_stack(reverse = T), color = "black", width = 0.7) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 8.5), breaks = 1:8) +
    scale_fill_manual(values = c("grey", rep("white", 20)) %>% setNames(1:21), breaks = c(1,2), labels = c("largest contig", "other contig")) +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "genome", y = "genome size (Mbp)", title = "all contigs")

p2 <- contigs_large %>%
    group_by(genome_id) %>%
    mutate(contig_ordered = factor(1:n())) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = contig_length/10^6, group = contig_ordered, fill = contig_ordered), position = position_stack(reverse = T), color = "black", width = 0.7) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 8.5), breaks = 1:8) +
    scale_fill_manual(values = c("grey", rep("white", 20)) %>% setNames(1:21), breaks = c(1,2), labels = c("largest contig", "other contig")) +
    theme_bw() +
    theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "genome", y = "genome size (Mbp)", title = "contigs > 500k")

p <- plot_grid(p1, p2, nrow = 2, align = "v", axis = "rl", labels = c("A", "B"), scale = 0.95) + paint_white_background()
ggsave(paste0(folder_data, "temp/12a-02-genome_size.png"), plot = p, width = 12, height = 8)

# genome size
g_size <- contigs_large %>%
    group_by(genome_id) %>%
    summarize(genome_size = sum(contig_length)/10^6)

mean(g_size$genome_size) # 6.93 Mbp

