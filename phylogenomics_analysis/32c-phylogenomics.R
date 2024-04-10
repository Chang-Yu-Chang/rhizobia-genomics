#' This scripts implements the phylogenetic comparative analysis

renv::load()
library(tidyverse)
library(janitor)
library(cowplot)
library(ape)
library(phytools)
library(geiger)
library(phangorn) # For rooting the tree
library(nlme) # For PGLS
source(here::here("analysis/00-metadata.R"))

tree_jaccard <- read.tree(paste0(folder_data, "temp/32-jaccard.tre"))
isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
isolates_contigs <- read_csv(paste0(folder_data, "temp/14-isolates_contigs.csv"))

# Use the 12 strains that has symbiosis data
isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
isolates_traits <- isolates_traits %>% drop_na(dry_weight_mg)

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
#cladelabels(tree_jaccard, c("clade 1","clade 2","clade 3"), c(47,41,29),wing.length=0,offset=0.5)
dev.off()

# Check names for 11 strains with both growth and symbiosis data
rownames(isolates_traits) <- isolates_traits$genome_id
chk <- name.check(tree_jaccard, isolates_traits)
tr <- drop.tip(tree_jaccard, chk$tree_not_data)
tb <- isolates_traits[tr$tip.label,]
rwn <- tb$genome_id
tb <- select(tb, starts_with("r_"), starts_with("lag_"), starts_with("maxOD_"), dry_weight_mg, root_weight_mg, nodule_number)
tb[is.na(tb)] <- 0
rownames(tb) <- rwn
name.check(tr, tb)

# PICs
# Compute PICs
t1 <- setNames(tb$r_30c, rownames(tb))
t2 <- setNames(tb$dry_weight_mg, rownames(tb))
pic_t1 <- pic(t1, tr)
pic_t2 <- pic(t2, tr)

# Fit a LM to PICs
lm_traits <- lm(dry_weight_mg ~ r_30c, data = tb) 
summary(lm_traits)
lm_pics <- lm(pic_t2 ~ pic_t1 + 0) 
summary(lm_pics)

# PGLS
gen <- rownames(tb)
cor_bm <- corBrownian(phy = tr, form = ~gen)
# for one trait
pgls_bm <- gls(dry_weight_mg ~ r_30c, data = tb, correlation = cor_bm)
summary(pgls_bm)

lm_pics$coefficients[1] # slope of PIC test -8.714665
pgls_bm$coefficients[2] # slope of PGLS -8.805389
abs(lm_pics$coefficients[1] -pgls_bm$coefficients[2])

tr
tb

# PGLS
cor_lambda <- corPagel(value = 1, phy = tr, form = ~gen)
# For one trait
pgls_lambda <- gls(dry_weight_mg ~ r_30c, data = tb, correlation = cor_lambda)
summary(pgls_lambda)
# For multiple traits
pgls_lambda <- gls(dry_weight_mg ~ r_30c + lag_30c + maxOD_30c, data = tb, correlation = cor_lambda)
anova(pgls_lambda)
pgls_lambda <- gls(dry_weight_mg ~ r_30c + r_25c + r_35c, data = tb, correlation = cor_lambda)
anova(pgls_lambda)



# Plot
## Orginial trait
p1 <- tb %>%
    ggplot() +
    geom_point(aes(x = r_30c, y = dry_weight_mg)) +
    geom_abline(intercept = lm_traits$coefficients[1], slope = lm_traits$coefficients[2]) +
    # geom_abline(intercept = 0, slope = lm_pics$coefficients[1], color = "maroon") +
    # geom_vline(xintercept = 0, linetype = 3) +
    # geom_hline(yintercept = 0, linetype = 3) +
    theme_classic()

## PICs
p2 <- tibble(pic1 = pic_t1, pic2 = pic_t2) %>%
    ggplot() +
    geom_point(aes(x = pic1, y = pic2)) +
    geom_abline(intercept = 0, slope = lm_pics$coefficients[1], color = "maroon") +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_hline(yintercept = 0, linetype = 3) +
    theme_classic()
p <- plot_grid(p1, p2, nrow = 1, scale = 0.9, align = "hv", axis = "tblr") + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "temp/32c-02-lm_vs_pic.png"), p, width = 6, height = 3)








# Use all 31 strains with growth data
# Check names for all strains
isolates_traits <- read_csv(paste0(folder_data, "temp/29-isolates_traits.csv"))
rownames(isolates_traits) <- isolates_traits$genome_id
chk <- name.check(tree_jaccard, isolates_traits)
tr <- drop.tip(tree_jaccard, chk$tree_not_data)
tb <- isolates_traits[tr$tip.label,]
rwn <- tb$genome_id
tb <- select(tb, starts_with("r_"), starts_with("lag_"), starts_with("maxOD_"))
tb[is.na(tb)] <- 0
rownames(tb) <- rwn
name.check(tr, tb)

# PGLS
rownames(tb)
gen <- rownames(tb)
cor_bm <- corBrownian(phy = tr, form = ~gen)
pgls_bm <- gls(lag_30c ~ r_30c, data = tb, correlation = cor_bm)
summary(pgls_bm)

# PGLS with lambda
cor_lambda <- corPagel(value = 1, phy = tr, form = ~gen)
pgls_lambda <- gls(lag_30c ~ r_30c, data = tb, correlation = cor_lambda)
summary(pgls_lambda)

# PGLS with multiple traits
pgls_ano <- gls(lag_30c ~ r_30c + r_25c, data = tb, correlation = cor_lambda)
anova(pgls_ano)
