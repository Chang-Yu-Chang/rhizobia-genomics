#' This script is to have a overview on the raw read statistics

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

compute_q <- function (asc) {
    #' Compute the mean phred (Quality score) of a raw read
    #' equation is from https://labs.epi2me.io/quality-scores/
    phred <- as.numeric(charToRaw(asc))-33
    round(-10*log10(sum(10^(-phred/10))/nchar(asc)), 2)
}

bash_var <- commandArgs(trailingOnly = TRUE)
file_reads = bash_var[1]
file_fig = bash_var[2]

#file_reads <- "/Users/cychang/Dropbox/lab/local-adaptation/data/genomics/genomes/Chang_Q5C_2/01-reads_qc/filtered_reads.txt"

# 0. read data ----
reads <- read_table(file_reads, col_names = c("name", "asc", "length"), show_col_types = F) %>%
    rowwise() %>%
    mutate(q_score = compute_q(asc), .keep = "unused") %>%
    ungroup()

median_q_score <- median(reads$q_score)

# 1. Read length vs. q score ----
p1 <- reads %>%
    ggplot() +
    geom_hline(yintercept = median_q_score, linetype = 2, color = "maroon") +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    geom_text(y = median_q_score, label = paste0("median:", median_q_score), x = Inf, color = "maroon", hjust = 1.2, vjust = -1) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(limits = c(0, 41)) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

# 2. histogram of read length ----
range(reads$length)
p2 <- reads %>%
    ggplot() +
    geom_histogram(aes(x = length/1000), color = "black", fill = "white", binwidth = 1) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)", y = "read count")


# 3. weighted histogram of read length ----
p3 <- reads %>%
    ggplot() +
    geom_histogram(aes(x = length/1000, weight = length/1000, y = after_stat(count/1000) ), color = "black", fill = "white", binwidth = 1) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)", y = " # of bases (kb)")

p <- plot_grid(p1,p2,p3, nrow = 2, align = "hv", scale = 0.96, labels = LETTERS[1:4]) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(filename = file_fig, plot = p, width = 10, height = 10)



















