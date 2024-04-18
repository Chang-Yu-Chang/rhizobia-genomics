#' This scripts implements the phylogenetic comparative analysis

renv::load()
library(tidyverse)
library(cowplot)
library(phytools)
library(phangorn) # For rooting the tree
library(ggtree)
library(tidytree)
library(geiger) # For checking names
#library(nlme) # For PGLS
source(here::here("metadata.R"))

# Traits
isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
isolates_traits <- read_csv(paste0(folder_data, "phenotypes_analysis/isolates_traits.csv"))
isolates <- isolates %>%
    left_join(isolates_contigs) %>%
    left_join(isolates_traits) %>%
    filter(!genome_id %in% c("g20", "g28"))

# Tree
load(file = paste0(folder_data, "phylogenomics_analysis/trees/trees.rdata"))
#tr, tr_acce, gpatl
# tr <- read.tree(paste0(folder_data, "genomics/mltree/isolates_core_b/aln.treefile"))
# list_others <- c(paste0("g", c(20, 28, 38:43)), "em1022", "usda1106", "em1021", "wsm419")
# tr <- tr %>% drop.tip(list_others)

## Compute phylogenetic signals
list_traits <- str_subset(names(isolates), "0c$|5c$|nodule|biomass|contig_length")
compute_ps <- function(tree, genome_id, trait_value) {
    set.seed(1)
    temp <- setNames(trait_value, genome_id)
    ps1 <- phylosig(tree, temp, test = T, nsim = 1000, method = "K")
    ps2 <- phylosig(tree, temp, test = T, nsim = 1000, method = "lambda")
    return(tibble(k = ps1$K, p_k = ps1$P, lambda = ps2$lambda, p_lambda = ps2$P))
}

tb1 <- tibble(trait = list_traits) %>%
    rowwise() %>%
    mutate(result = list(compute_ps(tr, isolates$genome_id, isolates[[trait]]))) %>%
    unnest(result)

tb2 <- tibble(trait = list_traits) %>%
    rowwise() %>%
    mutate(result = list(compute_ps(tr_acce, isolates$genome_id, isolates[[trait]]))) %>%
    unnest(result)

traits_ps <- bind_rows(mutate(tb1, tree = "core"), mutate(tb2, tree = "gcv"))
write_csv(traits_ps, paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))

# 1. Plot tree with traits ----
traits_ps <- read_csv(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))

p1 <- tr %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(color = species), align = T) +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top"
    )

ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)

p2 <- isolates %>%
    mutate(genome_id = factor(genome_id, ordered_tips)) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = r_30c)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2)
    ) +
    guides() +
    labs()

tmp <- traits_ps %>%
    mutate(across(c(k, p_k, lambda, p_lambda), function(x) round(x, 3))) %>%
    filter(trait == "r_30c", tree == "core")

p_ps <- tmp %>%
    ggplot() +
    geom_text(x = 0.5, y = 0.5, aes(label = paste("K = ", k, ", P = ", p_k, "\n\u03bb = ", lambda, ", P = ", p_lambda)), hjust = 0) +
    theme_classic()


p <- plot_grid(p1, p2, axis = "tb", align = "h")
ggsave(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/01-trait_ps.png"), p, width = 10, height = 6)

# 2. ----
traits_ps <- read_csv(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/traits_ps.csv"))

p1 <- tr_acce %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(color = species), align = T) +
    scale_color_manual(values = species_colors) +
    theme_tree() +
    theme(
        legend.position = "top"
    )

ordered_tips <- p1$data %>%
    filter(isTip) %>%
    arrange(y) %>%
    pull(label)

p2 <- isolates %>%
    mutate(genome_id = factor(genome_id, ordered_tips)) %>%
    ggplot() +
    geom_col(aes(x = genome_id, y = r_30c)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2)
    ) +
    guides() +
    labs()

tmp <- traits_ps %>%
    mutate(across(c(k, p_k, lambda, p_lambda), function(x) round(x, 3))) %>%
    filter(trait == "r_30c", tree == "core")

p_ps <- tmp %>%
    ggplot() +
    geom_text(x = 0.5, y = 0.5, aes(label = paste("K = ", k, ", P = ", p_k, "\n\u03bb = ", lambda, ", P = ", p_lambda)), hjust = 0) +
    theme_classic()


p <- plot_grid(p1, p2, axis = "tb", align = "h")
ggsave(paste0(folder_data, "phylogenomics_analysis/phylogenetic_signals/01-trait_ps.png"), p, width = 10, height = 6)




if (FALSE) {

# Root the tree
class(tr)
tr <- midpoint(tr)

names(tr)
tr$tip.label
plot(tr)
tb <- isolates %>%
    mutate(genome_id = factor(genome_id,  tr$tip.label)) %>%
    arrange(genome_id) %>%
    left_join(distinct(isolates_contigs, genome_id, species) %>% select(genome_id, species))


# Plot
tr %>%
    as_tibble() %>%
    left_join(rename(isolates, label = genome_id)) %>%
    as.treedata() %>%
    ggtree() +
    geom_tiplab(aes(label = label))


# Check names for 11 strains with both growth and symbiosis data
rownames(isolates) <- isolates$genome_id
chk <- name.check(tr, isolates)
tr <- drop.tip(tr, chk$tree_not_data)
tb <- isolates[tr$tip.label,]
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
isolates <- read_csv(paste0(folder_data, "temp/29-isolates.csv"))
rownames(isolates) <- isolates$genome_id
chk <- name.check(tr, isolates)
tr <- drop.tip(tr, chk$tree_not_data)
tb <- isolates[tr$tip.label,]
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
}
