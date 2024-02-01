#' This script plot the NJ trees

renv::load()
library(tidyverse)
library(cowplot)
library(ggsci)
library(ape)
library(ggtree)
library(tidytree)
source(here::here("analysis/00-metadata.R"))

load(paste0(folder_data, "temp/32-trees.RData"))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))
isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
isolates_contigs <- isolates_contigs %>%
    filter(!genome_id %in% c("g20", "g28", "g2", "g3", "g15"))  
list_traits <- c("dry_weight_mg", "nodule_number", "root_weight_mg",
        str_subset(colnames(isolates_traits), "r_"),
        str_subset(colnames(isolates_traits), "lag_"),
        str_subset(colnames(isolates_traits), "maxOD_"))
dists <- read_csv(paste0(folder_data, 'temp/31-dists.csv'))

# dist_genetics <- read_csv(paste0(folder_data, "temp/19-dist_genetics.csv"))
# dist_traits <- read_csv(paste0(folder_data, "temp/29-dist_traits.csv"))
# dists_long <- read_csv(paste0(folder_data, "temp/31-dists_long.csv"))



# 1. Plot the trees by populations
plot_tree1 <- function (tree, gtitle) {
    tree %>%
        as_tibble() %>%
        left_join(rename(isolates_traits, label = genome_id), by = "label") %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(aes(color = population), size = 3) +
        #geom_tippoint(aes(color = population)) +
        scale_color_manual(values = c(MLBS = "#0F6290", Phila = "#85325C")) +
        theme_tree(legend.position = 'centre') +
        theme(
            legend.position = c(0.85,0.25),
            legend.background = element_rect(fill = NA, color = NA)
        ) +
        #guides(color = guide_legend(override.aes = aes(label = ""), title = NULL)) +
        guides(color = guide_legend(title = NULL)) +
        labs(title = gtitle)
}
p_trees <- rep(list(NA), length(list_trees))
names(p_trees) <- names(list_trees)
for (i in 1:length(list_trees)) p_trees[[i]] <- plot_tree1(list_trees[[i]], names(list_trees)[i])
p <- plot_grid(plotlist = p_trees, nrow = 2, scale = 0.85) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-01-trees_population.png"), p, width = 20, height = 10)


# 2. Plot the trees by species
plot_tree2 <- function (tree, gtitle) {
    tree %>%
        as_tibble() %>%
        left_join(rename(isolates_contigs, label = genome_id), by = "label") %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(aes(color = species), size = 3) +
        #geom_tippoint(aes(color = species)) +
        scale_color_manual(values = c(meliloti = "#423E37", medicae = "#E3B23C")) +
        theme_tree(legend.position = 'centre') +
        theme(
            legend.position = c(0.85,0.25),
            legend.background = element_rect(fill = NA, color = NA)
        ) +
        #guides(color = guide_legend(override.aes = aes(label = ""), title = NULL)) +
        guides(color = guide_legend(title = NULL)) +
        labs(title = gtitle)
}
p_trees <- rep(list(NA), length(list_trees))
names(p_trees) <- names(list_trees)
for (i in 1:length(list_trees)) p_trees[[i]] <- plot_tree2(list_trees[[i]], names(list_trees)[i])
p <- plot_grid(plotlist = p_trees, nrow = 2, scale = 0.85) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-02-trees_species.png"), p, width = 20, height = 10)


# 3. Plot trees for each traits, colored by population
make_tree <- function(di, d_trait) {
    #' This functions creates a phylo object from a long-format distance matrix
    if (d_trait %in% paste0("d_", list_traits)) {
        list_avails <- isolates_traits$genome_id[!is.na(isolates_traits[str_remove(d_trait, "d_")])]
        di <- bind_cols(select(di, genome_id1, genome_id2), tibble(dd = unlist(di[,d_trait]))) %>%
            filter(genome_id1 %in% list_avails, genome_id2 %in% list_avails)
    } else {
        di <- bind_cols(select(di, genome_id1, genome_id2), tibble(dd = unlist(di[,d_trait]))) 
    }
    colnames(di)[3] <- "dd"
    dist1 <- tibble(genome_id1 = di$genome_id1, genome_id2 = di$genome_id2, dd = di$dd) 
    dist1_swapped <- tibble(genome_id1 = di$genome_id2, genome_id2 = di$genome_id1, dd = di$dd)
    dist1_swapped <- filter(dist1_swapped, genome_id1 != genome_id2)
    dist2 <- bind_rows(dist1, dist1_swapped) %>%
        drop_na(dd) %>%
        pivot_wider(names_from = genome_id2, values_from = dd) %>%
        select(-genome_id1) %>%
        as.matrix() 
    tree <- dist2 %>%
        nj() 
    return(tree)
}
p_trees <- rep(list(NA), length(list_traits))
for(i in 1:length(list_traits)) p_trees[[i]] <- dists %>% make_tree(paste0("d_", list_traits[i])) %>% plot_tree1(list_traits[i])
p <- plot_grid(plotlist = p_trees, nrow = 3, scale = 0.9) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-03-trees_trait_population.png"), p, width = 20, height = 12)

# 4. Plot trees for each traits, colored by species
p_tree <- rep(list(NA), length(list_traits))
for(i in 1:length(list_traits)) p_tree[[i]] <- dists %>% make_tree(paste0("d_", list_traits[i])) %>% plot_tree2(list_traits[i])
p <- plot_grid(plotlist = p_tree, nrow = 3, scale = 0.9) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-04-trees_trait_species.png"), p, width = 20, height = 12)

# 5. Plot the meliloti trees
plot_tree3 <- function (tree, gtitle) {
    tree %>%
        as_tibble() %>%
        left_join(rename(isolates_contigs, label = genome_id), by = "label") %>%
        as.treedata() %>%
        ggtree() +
        geom_tiplab(size = 3) +
        #geom_tippoint(aes(color = species)) +
        #scale_color_manual(values = c(meliloti = "#423E37", medicae = "#E3B23C")) +
        theme_tree(legend.position = 'centre') +
        theme(
            legend.position = c(0.85,0.25),
            legend.background = element_rect(fill = NA, color = NA)
        ) +
        #guides(color = guide_legend(override.aes = aes(label = ""), title = NULL)) +
        guides(color = guide_legend(title = NULL)) +
        labs(title = gtitle)
}

p_trees <- rep(list(NA), length(list_trees_meliloti))
names(p_trees) <- names(list_trees_meliloti)
for (i in 1:length(list_trees_meliloti)) p_trees[[i]] <- plot_tree3(list_trees_meliloti[[i]], names(list_trees_meliloti)[i])
p <- plot_grid(plotlist = p_trees, nrow = 2, scale = 0.85) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-05-trees_meliloti.png"), p, width = 20, height = 10)

# 6. Plot the medicae trees
p_trees <- rep(list(NA), length(list_trees_medicae))
names(p_trees) <- names(list_trees_meliloti)
for (i in 1:length(list_trees_medicae)) p_trees[[i]] <- plot_tree3(list_trees_medicae[[i]], names(list_trees_medicae)[i])
p <- plot_grid(plotlist = p_trees, nrow = 2, scale = 0.85) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32a-06-trees_medicae.png"), p, width = 20, height = 10)
