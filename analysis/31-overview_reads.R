#' This script is to have a overview on the raw read statistics

library(tidyverse)
library(janitor)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. read the raw read data ----
if (FALSE) {
    #' Aggregate the raw read data
    #' The trunck only needs to run once
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
    mutate(genome_id = factor(paste0("g", genome_id), paste0("g", 1:19)))

write_csv(raw_reads, file = paste0(folder_data, "temp/31-raw_reads.csv"))

}
raw_reads <- read_csv(paste0(folder_data, "temp/31-raw_reads.csv"), show_col_types = F)
isolates_rhizo <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F)

# Clean up
raw_reads <- raw_reads %>% mutate(genome_id = factor(genome_id, paste0("g", 1:19)))

# Summary stat for each file number histogram
raw_reads_summ <- raw_reads %>%
    group_by(genome_id) %>%
    summarize(n = n(),
              median_q_score = median(q_score),
              median_length = median(length))

# Remove rhizobium and keep only Ensifer
raw_reads <- raw_reads %>%
    filter(genome_id %in% isolates_rhizo$genome_id[isolates_rhizo$genus == "Ensifer"])

# 1. genome 1 read length vs. q_score ----
p <- raw_reads %>%
    filter(genome_id == "g2") %>%
    ggplot() +
    geom_hline(data = raw_reads_summ[1,], aes(yintercept = median_q_score), linetype = 2, color = "maroon") +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(limits = c(0, 41)) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")
ggsave(paste0(folder_data, "temp/31-01-length_vs_q.png"), plot = p, width = 4, height = 4)

# 2. read length vs. q_score, facet ----
p <- raw_reads %>%
    ggplot() +
    geom_hline(data = raw_reads_summ, aes(yintercept = median_q_score), linetype = 2, color = "maroon") +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    geom_text(data = raw_reads_summ, aes(y = median_q_score, label = paste0("median:", median_q_score)), x = Inf, color = "maroon", hjust = 1.2, vjust = -1) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(limits = c(0, 41)) +
    facet_wrap(~genome_id, ncol = 4, scales = "free_x") +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "temp/31-02-length_vs_q_facet.png"), plot = p, width = 10, height = 10)

# 3. read length vs. q_score, density plot ----
p <- raw_reads %>%
    filter(genome_id == "g2") %>%
    ggplot() +
    geom_density2d_filled(aes(x = length/1000, y = q_score), alpha = 0.5) +
    scale_x_continuous(labels = scales::label_comma(), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 41)) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "temp/31-03-length_vs_q_density.png"), plot = p, width = 6, height = 6)

# 4. read length vs. q_score, density plot and facet ----
p <- raw_reads %>%
    ggplot() +
    geom_density2d_filled(aes(x = length/1000, y = q_score), alpha = 0.5) +
    scale_x_continuous(labels = scales::label_comma(), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 41)) +
    facet_wrap(~genome_id, ncol = 3, nrow = 5, scales = "free_x") +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "temp/31-04-length_vs_q_density_facet.png"), plot = p, width = 10, height = 15)

# 5. plot read length for all reads. violin ----
p <- raw_reads %>%
    ggplot() +
    geom_violin(aes(x = genome_id, y = length/1000), draw_quantiles = 0.5) +
    scale_y_continuous(labels = scales::label_comma(), trans = "log10", breaks = 10^(c(-2,-1,0,1,2,3))) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = "grey", linetype = 2)
    ) +
    guides() +
    labs(x = "sample", y = "read length (kb)")

ggsave(paste0(folder_data, "temp/31-05-read_length.png"), plot = p, width = 6, height = 4)

# 6. plot reads below 80k ----
raw_reads_80k_summ <- raw_reads %>%
    filter(length < 80000) %>%
    group_by(genome_id) %>%
    summarize(n = n(),
              median_q_score = median(q_score),
              median_length = median(length))
p <- raw_reads %>%
    filter(length < 80000) %>%
    ggplot() +
    geom_hline(data = raw_reads_summ, aes(yintercept = median_q_score), linetype = 2, color = "maroon") +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    geom_text(data = raw_reads_80k_summ, aes(y = median_q_score, label = paste0("median:", median_q_score)), x = Inf, color = "maroon", hjust = 1.2, vjust = -1) +
    scale_x_continuous(labels = scales::label_comma(), expand = c(0,0), limits = c(0, 80)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 41)) +
    facet_wrap(~genome_id, ncol = 4, nrow = 4) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "temp/31-06-length_vs_q_facet_80k.png"), plot = p, width = 10, height = 10)























