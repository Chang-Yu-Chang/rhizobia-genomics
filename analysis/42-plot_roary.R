#' This script plot the pangenome analysis

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 0. read data
gpa_meta <- read_csv("~/mgjw/problem_set7/roary/gene_presence_absence.csv", show_col_types = F)
gpa <- read_table("~/mgjw/problem_set7/roary/gene_presence_absence.Rtab", show_col_types = F) %>% select(-`NC_017512.1:`)

# pivot longer
gpa[gpa == 0] <- NA
gpa_long <- gpa %>%
    pivot_longer(cols = -Gene, values_drop_na = T) %>%
    clean_names()

# 1. some numbers
# Size of the Core Genome

# Size of the Accessory Genome
# Pan Genome Size


# 2. Number of gene shared in these isolates
gpa_long %>%
    group_by(gene) %>%
    count(name = "n_isolates") %>%
    mutate(n_isolates = factor(n_isolates, 1:8)) %>%
    group_by(n_isolates) %>%
    count(name = "n_genes") %>%
    ggplot() +
    geom_col(aes(x = n_isolates, y = n_genes), color = "black", fill = "white") +
    #scale_x_continuous(breaks = 1:8) +
    theme_classic() +
    theme() +
    guides() +
    labs()

# 3. heatmap for genes
gpa_long %>%
    filter(str_detect(gene, "b")) %>%
    mutate(value = factor(value)) %>%
    ggplot() +
    geom_tile(aes(x = name, y = gene, fill = value)) +
    scale_fill_manual(values = c(`1` = "maroon")) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "strain", y = "gene")




















