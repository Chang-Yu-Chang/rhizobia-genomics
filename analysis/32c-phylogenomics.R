#' This scripts implements the phylogenetic comparative analysis

renv::load()
library(tidyverse)
library(janitor)
library(ape)
#library(tidytree)
library(phytools)
library(geiger)
library(phangorn)
source(here::here("analysis/00-metadata.R"))

tree_jaccard <- read.tree(paste0(folder_data, "temp/32-jaccard.tre"))
isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))

# Root the tree
class(tree_jaccard)
tree_jaccard <- midpoint(tree_jaccard)

names(tree_jaccard)
tree_jaccard$tip.label

tb <- isolates_traits %>%
    filter(genome_id %in% tree_jaccard$tip.label) %>%
    mutate(genome_id = factor(genome_id,  tree_jaccard$tip.label)) %>%
    arrange(genome_id) %>%
    left_join(distinct(isolates_contigs, genome_id, species) %>% select(genome_id, species))

# Plot
png(filename = paste0(folder_data, "temp/32c-01-tree_jaccard.png"), width = 20, height = 20, units = "cm", res = 1000)
plotTree(tree_jaccard, fsize = 1, lwd = 2, split.vertical = TRUE)
nodelabels(bg = "white", cex = 0.5, frame = "circle")
tiplabels(tb$genome_id, col = c("blue", "red")[match(tb$species, c("meliloti", "medicae"))], 
    bg = NA, frame = NULL, offset = 5)
#cladelabels(tree_jaccard, c("clade 1","clade 2","clade 3"), c(47,41,29),wing.length=0,offset=0.5)
dev.off()

# Check names
rownames(isolates_traits) <- isolates_traits$genome_id
chk <- name.check(tree_jaccard, isolates_traits)
tr <- drop.tip(tree_jaccard, chk$tree_not_data)
tb <- isolates_traits[tr$tip.label,]
rwn <- tb$genome_id
tb <- select(tb, dry_weight_mg, root_weight_mg, nodule_number)
#tb[is.na(tb)] <- 0
rownames(tb) <- rwn
name.check(tr, tb)

# Phylo PCA
pca <- phyl.pca(tr, tb)
scores(pca)
phylomorphospace(tr, scores(pca)[,1:2])

# 
tb %>%
    ggplot() +
    geom_point(aes(x = r_30c, y = r_25c)) +
    theme_classic()

# Compute PICs
t1 <- setNames(tb$r_30c, rownames(tb))
t2 <- setNames(tb$r_25c, rownames(tb))
pic_t1 <- pic(t1, tr)
pic_t2 <- pic(t2, tr)

# Fit a LM to contrasts
lm(r_25c ~ r_30c, data = tb) %>% summary()
lm(pic_t2 ~ pic_t1 + 0) %>% summary()
