#' This script plot the Fst per gene

renv::load()
library(tidyverse)
library(ggsci)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
set1_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/set1_fst.csv"))
set1_pop_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/set1_pop_fst.csv"))


# 1. Plot the Fst within gradients ----
n_sccg <- distinct(set1_fst, gene, singlecopy) %>% pull(singlecopy) %>% sum
p <- set1_fst %>%
    filter(metric == "Gst_est", singlecopy) %>%
    select(-metric, -singlecopy) %>%
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
    labs(x = "Nei's G_st", title = paste0(n_sccg, " single-copy core genes"))
ggsave(paste0(folder_data, "genomics_analysis/fst/01-fst_sccg.png"), p, width = 4, height = 4)


# 2. Plot the Fst between all populations across the two gradients ----
n_sccg <- distinct(set1_fst, gene, singlecopy) %>% pull(singlecopy) %>% sum
p <- set1_pop_fst %>%
    filter(singlecopy) %>%
    select(-singlecopy) %>%
    rename(gradient1 = pop1, gradient2 = pop2) %>%
    #pivot_longer(cols = -c(gene, Gst_est), names_to = "gradient_group", values_to = "pop_group") %>%
    ggplot() +
    geom_histogram(aes(x = Gst_est), alpha = .5, position = "identity", binwidth = 0.01) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    scale_fill_aaas() +
    coord_cartesian(clip = "off") +
    facet_grid(gradient1 ~ gradient2) +
    theme_bw() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(.8,.8),
        legend.background = element_rect(color = "grey10", fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
    ) +
    guides() +
    labs(x = "Nei's G_st", title = paste0(n_sccg, " single-copy core genes"))
ggsave(paste0(folder_data, "genomics_analysis/fst/02-fst_pop_sccg.png"), p, width = 8, height = 8)


# Top genes
set1_fst %>%
    filter(metric == "Gst_est", singlecopy) %>%
    select(-metric, -singlecopy) %>%
    pivot_longer(-gene, names_to = "gradient") %>%
    mutate(gradient = ifelse(gradient == "elev", "elevation", "urbanization")) %>%
    group_by(gradient) %>%
    arrange(desc(value)) %>%
    slice(1:ceiling(n()/100))











