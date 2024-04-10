#' This script summarize assembled genome information

renv::load()
library(tidyverse)
library(janitor)
library(seqinr)
source(here::here("metadata.R"))

genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
filtered_reads <- read_csv(paste0(folder_data, "genomics_analysis/raw_reads/filtered_reads.csv"))

# Remove reads that are too long
sum(filtered_reads$length > 80000) # 9 reads with > 80kb length
filtered_reads <- filtered_reads %>%
    mutate(genome_id = factor(genome_id, genomes$genome_id)) %>%
    filter(length < 80000)

# Number of reads per sample
n_reads <- filtered_reads %>%
    group_by(genome_id) %>%
    dplyr::count(name = "n_reads") %>%
    ungroup()
range(n_reads$n_reads) # 47161 234638
median(n_reads$n_reads) # 104054.5


# 1. Read length vs. q score
filtered_reads_median <- filtered_reads %>%
    group_by(genome_id) %>%
    summarize(median_q_score = median(q_score), median_length = median(length))

p <- filtered_reads %>%
    ggplot() +
    geom_hline(data = filtered_reads_median, aes(yintercept = median_q_score), linetype = 2, color = "maroon") +
    geom_point(aes(x = length/1000, y = q_score), size = 0.2, alpha = 0.2) +
    geom_text(data = filtered_reads_median, aes(y = median_q_score, label = paste0("median:", median_q_score)), x = Inf, color = "maroon", hjust = 1.2, vjust = -1) +
    facet_wrap(~genome_id, ncol = 5) +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 1)
    ) +
    guides() +
    labs(x = "read length (kb)", y = "mean Phred")

ggsave(paste0(folder_data, "genomics_analysis/raw_reads/01-length_vs_q.png"), p, width = 10, height = 10)

# 2. histogram of read length
p <- filtered_reads %>%
    ggplot() +
    geom_histogram(aes(x = length/1000), color = "black", fill = "white", binwidth = 1) +
    facet_wrap(~genome_id, ncol = 5) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)", y = "read count")

ggsave(paste0(folder_data, "genomics_analysis/raw_reads/02-read_length.png"), p, width = 10, height = 10)

#
range(filtered_reads_median$median_length) # 2748 5367

# 3. weighted histogram of read length
p <- filtered_reads %>%
    ggplot() +
    geom_histogram(aes(x = length/1000, weight = length/1000, y = after_stat(count/1000) ), color = "black", fill = "white", binwidth = 1) +
    facet_wrap(~genome_id, ncol = 5) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "read length (kb)", y = " # of bases (kb)")

ggsave(paste0(folder_data, "genomics_analysis/raw_reads/03-weighted_read_length.png"), p, width = 10, height = 10)

#
range(filtered_reads_median$median_q_score) # 18.26 20.24

# 4. Calculat estimated coverage
filtered_reads_coverage <- filtered_reads %>%
    group_by(genome_id) %>%
    summarize(total_length = sum(length)) %>%
    mutate(coverage = total_length/(7*10^6))

round(range(filtered_reads_coverage$coverage)) # 48 171








