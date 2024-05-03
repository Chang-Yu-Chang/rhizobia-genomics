#' This script plots the Robinson-Foulds distance

renv::load()
library(tidyverse)
library(ape)
source(here::here("metadata.R"))

load(paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))

# Normalized Robinson-Foulds tree distance
dist.topo(tr, drop.tip(tr_gpa, "g20"), method = "PH85") %>% suppressWarnings()

# Random trees
trs <- rmtopology(1000, 30, tip.label = tr$tip.label)
trs <- c(trs, tr, drop.tip(tr_gpa, "g20"))
md <- dist.topo(trs)
dim(as.matrix(md)) # 1002

tbb <- as_tibble(as.matrix(md)/max(md)[1]) %>%
    mutate(row = paste0("tree", 1:1002)) %>%
    pivot_longer(cols = -row, names_to = "col", values_drop_na = T)

obs <- tbb %>% filter(row == "tree1001", col == "tree1002") %>% pull(value)

p <- tbb %>%
    #filter(row != col) %>%
    filter(row < col) %>%
    filter(row != "tree1001") %>%
    arrange(value) %>%
    ggplot() +
    geom_histogram(aes(x = value), breaks = seq(0, 1, 0.05), color = "black", fill = "white") +
    geom_vline(xintercept = obs, color = "maroon", linetype = 2) +
    #annotate("text", label = paste0("between core gene and\ngene content trees\n", round(obs,2)), x = 0.5, y = Inf, vjust = 1, hjust = 0.5, color = "maroon") +
    scale_y_log10() +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x =  "Normalized Robinson-Foulds Distance", y = "Count")

ggsave(here::here("plots/FigS9.png"), p, width = 4, height = 3)


#dist.topo(tr, rtree(30, rooted = T, tip.label = tr$tip.label)) %>% suppressWarnings()
#dist.topo(tr, drop.tip(tr_gpa, "g20"), method = "score") %>% suppressWarnings()
#dist.topo(tr, rtree(30, rooted = T, tip.label = tr$tip.label), method = "score") %>% suppressWarnings()
