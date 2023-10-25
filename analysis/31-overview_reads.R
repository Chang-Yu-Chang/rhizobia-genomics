#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data
if (FALSE) {
list_g <- paste0("Chang_Q5C_", 1:19)
list_reads <- rep(list(NA), 19)
compute_q <- function (asc) {
    #' Compute the mean phred (Quality score) of a raw read
    #' equation is from https://labs.epi2me.io/quality-scores/
    phred <- as.numeric(charToRaw(asc))-33
    round(-10*log10(sum(10^(-phred/10))/nchar(asc)), 2)
}

for (i in 1:19) {
    cat("\n", i)
    list_reads[[i]] <- read_table(paste0(folder_data, "temp/plasmidsaurus/", list_g[i], "/01-filtlong/raw_reads.txt"), col_names = c("name", "asc", "length"), show_col_types = F) %>%
        rowwise() %>%
        mutate(q_score = compute_q(asc), .keep = "unused")
}

raw_reads <- list_reads %>%
    lapply(ungroup) %>%
    bind_rows(.id = "genome_id") %>%
    mutate(genome_id = factor(genome_id, 1:19))

write_csv(raw_reads, file = paste0(folder_data, "temp/31-raw_reads.csv"))

}

raw_reads <- read_csv(paste0(folder_data, "temp/31-raw_reads.csv"))

# Summary stat for each file number histogram
raw_reads %>%
    group_by(genome_id) %>%
    summarize(n = n(),
              median_q_score = median(q_score),
              median_length = median(length))


# 1. genome 1 read length vs. q_score
p <- raw_reads %>%
    filter(genome_id == 1) %>%
    ggplot() +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    scale_x_continuous(labels = scales::label_comma()) +
    #scale_x_log10() +
    #facet_wrap(~genome_id, ncol = 5, nrow = 4) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")
ggsave(paste0(folder_data, "temp/31-01-g1_read_length_vs_q_score.png"), plot = p, width = 4, height = 4)

# 2. read length vs. q_score, facet
p <- raw_reads %>%
    ggplot() +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    scale_x_continuous(labels = scales::label_comma()) +
    facet_wrap(~genome_id, ncol = 5, nrow = 4, scales = "free_x") +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "temp/31-02-read_length_vs_q_score.png"), plot = p, width = 16, height = 20)

# 2a. read length vs. q_score, overlapping
p <- raw_reads %>%
    ggplot() +
    geom_point(aes(x = length/1000, y = q_score, color = genome_id), size = 0.2, alpha = 0.2) +
    scale_x_continuous(labels = scales::label_comma()) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "temp/31-02a-read_length_vs_q_score.png"), plot = p, width = 6, height = 6)


# 3. plot read length for all reads
p <- raw_reads %>%
    filter(genome_id==1) %>%
    ggplot() +
    geom_histogram(aes(x = length/1000), color = "black", fill = "white") +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)")

ggsave(paste0(folder_data, "temp/31-03-read_length.png"), plot = p, width = 4, height = 4)

# 4. genome 1, plot read length for reads below 20k
p <- raw_reads %>%
    filter(length < 20000) %>%
    filter(genome_id==1) %>%
    ggplot() +
    geom_histogram(aes(x = length/1000), color = "black", fill = "white") +
    scale_y_continuous(labels = scales::label_comma()) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)")

ggsave(paste0(folder_data, "temp/31-04-g1_read_length_20k.png"), plot = p, width = 4, height = 4)


# 5. plot read length for reads below 20k, facet
p <- raw_reads %>%
    filter(length < 20000) %>%
    ggplot() +
    geom_histogram(aes(x = length/1000), color = "black", fill = "white") +
    scale_y_continuous(labels = scales::label_comma()) +
    facet_wrap(~genome_id, ncol = 4, scales = "free_y") +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)")

ggsave(paste0(folder_data, "temp/31-05-read_length_20k.png"), plot = p, width = 16, height = 20)


# 5a. plot read length for reads below 20k, overlap
p <- raw_reads %>%
    filter(length < 20000) %>%
    ggplot() +
    geom_histogram(aes(x = length/1000, fill = genome_id), color = "black", position = position_identity(), alpha = 0.1) +
    scale_y_continuous(labels = scales::label_comma()) +
    #facet_wrap(~genome_id, ncol = 4, scales = "free_y") +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)")

ggsave(paste0(folder_data, "temp/31-05a-read_length_20k.png"), plot = p, width = 6, height = 6)






















