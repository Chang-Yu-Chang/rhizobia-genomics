#' This script generates figure 3

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
library(vegan) # for computing jaccard
source(here::here("analysis/00-metadata.R"))

site_group_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", 
                 `suburban` = "#0cc45f", `urban` = "#a642bf", control = "grey")


# Read gene presence-absence data
gpa <- read_csv(paste0(folder_data, "temp/13-gpa.csv"))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))

gpa %>%
    rename(genome_id = name) %>%
    pivot_longer(-genome_id) %>%
    group_by(name) %>%
    summarize(n_genomes = sum(value)) %>%
    filter(n_genomes < 41)


# Panel A. Cartoons
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png")) + draw_text("placeholder")



# Panel B. SNPs
snp <- read_table(paste0(folder_data, "genomics/variants/em1021_snippy/core.tab"))
gs <- read_table(paste0(folder_data, "genomics/variants/em1021_snippy/core.txt"))
vcf <- read.vcfR(paste0(folder_data, "genomics/variants/em1021_snippy/core.vcf"))

vcfR::getCHROM(vcf) %>% table() # SNPs on chromosome
gl <- vcfR2genlight(vcf) # VCF to genlight
snps <- tab(gl); dim(snps)


## Test VA populaitons: high vs low elevation 
isolates_i1 <- isolates_test %>% filter(population == "VA")
m <- snps[match(isolates_i1$genome_id, gl@ind.names),]; dim(m) # Filter the snps table for focal genomes
#table(apply(m, 2, sum))
dist_m1 <- vegdist(m, "jaccard")
permanova1 <- adonis2(dist_m1 ~ site_group, data = isolates_i1)
permanova1 # no
pcoa1 <- cmdscale(dist_m1, eig = T)
eigs1 <- round(pcoa1$eig / sum(pcoa1$eig)*100, 2)

## Test PA populaitons: high vs low elevation 
isolates_i2 <- isolates_test %>% filter(population == "PA")
m <- snps[match(isolates_i2$genome_id, gl@ind.names),]; dim(m) # Filter the snps table for focal genomes
#m <- m[,apply(m, 2, sum) != 0]; dim(m) # 323 snps??
dist_m2 <- vegdist(m, "jaccard")
permanova2 <- adonis2(dist_m2 ~ site_group, data = isolates_i2)
permanova2 # no
pcoa2 <- cmdscale(dist_m2, eig = T)
eigs2 <- round(pcoa2$eig / sum(pcoa2$eig)*100, 2)

## Plot the mds
p2_1 <- isolates_i1 %>% 
    bind_cols(tibble(mds1 = pcoa1$points[,1], mds2 = pcoa1$points[,2])) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group), size = 3, shape = 21, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.margin = unit(c(0,0,0,0), "mm"),
        #legend.background = element_rect(color = "black", fill = "white"),
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs1[1], "%)"), y = paste0("PCoA axis 1(", eigs1[2], "%)"))

p2_2 <- isolates_i2 %>% 
    bind_cols(tibble(mds1 = pcoa2$points[,1], mds2 = pcoa2$points[,2])) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group), size = 3, shape = 21, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.margin = unit(c(0,0,0,0), "mm"),
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs2[1], "%)"), y = paste0("PCoA axis 1(", eigs2[2], "%)"))







# Panel C. GCVs
isolates_test <- isolates %>% filter(genome_id %in% gpa$name, !is.na(site_group)) 
gpa_test <- gpa %>% filter(name %in% isolates_test$genome_id)

## Test VA populaitons: high vs low elevation 
isolates_i1 <- isolates_test %>% filter(population == "VA")
labels <- isolates_i1$site_group; labels <- as.factor(labels)
gpa_i <- gpa_test %>% filter(name %in% isolates_i1$genome_id)
m <- as.matrix(gpa_i[,-1]); dim(m)
m <- m[,apply(m, 2, sum) != 0]; dim(m) # 23575 genes
dist_m1 <- vegdist(m, "jaccard")
permanova1 <- adonis2(dist_m1 ~ site_group, data = isolates_i1)
permanova1 # no
pcoa1 <- cmdscale(dist_m1, eig = T)
eigs1 <- round(pcoa1$eig / sum(pcoa1$eig)*100, 2)

## Test PA populaitons: high vs low elevation 
isolates_i2 <- isolates_test %>% filter(population == "PA")
labels <- isolates_i2$site_group; labels <- as.factor(labels)
gpa_i <- gpa_test %>% filter(name %in% isolates_i2$genome_id)
m <- as.matrix(gpa_i[,-1]); dim(m)
m <- m[,apply(m, 2, sum) != 0]; dim(m) # 12946 genes
dist_m2 <- vegdist(m, "jaccard")
permanova2 <- adonis2(dist_m2 ~ site_group, data = isolates_i2)
permanova2 # no
pcoa2 <- cmdscale(dist_m2, eig = T)
eigs2 <- round(pcoa2$eig / sum(pcoa2$eig)*100, 2)

## Plot the mds
p3_1 <- isolates_i1 %>% 
    bind_cols(tibble(mds1 = pcoa1$points[,1], mds2 = pcoa1$points[,2])) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group), size = 3, shape = 21, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",
        # legend.margin = unit(c(0,0,0,0), "mm"),
        # legend.background = element_rect(color = "black", fill = "white"),
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs1[1], "%)"), y = paste0("PCoA axis 1(", eigs1[2], "%)"))

p3_2 <- isolates_i2 %>% 
    bind_cols(tibble(mds1 = pcoa2$points[,1], mds2 = pcoa2$points[,2])) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey80") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey80") +
    geom_point(aes(x = mds1, y = mds2, color = site_group), size = 3, shape = 21, stroke = 1) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none",
        #legend.margin = unit(c(0,0,0,0), "mm"),
        #legend.position = c(0.25,0.25),
        #legend.background = element_rect(color = "black", fill = "white"),
        legend.title = element_blank()
    ) +
    guides() +
    labs(x = paste0("PCoA axis 1(", eigs2[1], "%)"), y = paste0("PCoA axis 1(", eigs2[2], "%)"))


# 
p_right <- plot_grid(p2_1, p2_2, p3_1, p3_2, nrow = 2, rel_heights = c(1, 1), scale = 0.95,
    align = "vh", axis = "tblr", labels = c("B", "", "C", ""))

p <- plot_grid(p1, p_right, rel_widths = c(1,1.5), labels = c("A", "")) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 6)
 




if (FALSE) {
library(randomForest)

isolates_test <- isolates %>% filter(genome_id %in% gpa$name, !is.na(site_group)) 
gpa_test <- gpa %>% filter(name %in% isolates_test$genome_id)

# Test VA populaitons: high vs low elevation 
isolates_i <- isolates_test %>% filter(population == "VA")
labels <- isolates_i$site_group; labels <- as.factor(labels)
gpa_i <- gpa_test %>% filter(name %in% isolates_i$genome_id)
m <- as.matrix(gpa_i[,-1]); dim(m)
rf_model <- randomForest(m, y = labels)
print(paste("Out-of-bag (OOB) error rate:", rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]))
print(table(labels, rf_model$predicted))
print(sort(importance(rf_model), decreasing = T)[1:10])




gpa %>% 
    pca <- glPca(gl, nf = 10, center = T, scale = T)
# convert scores of vcf.pca into a tibble
sc <- as_tibble(pca$scores) %>%
    bind_cols(isolates_arr)
pcs <- round(100 * pca$eig / sum(pca$eig), 2)
p <- tibble(eigenvalue = 1:length(pca$eig), pvar = pcs) %>%
    ggplot() +
    geom_col(aes(x = eigenvalue, y = pvar)) +
    scale_x_continuous(limits = c(0.5,10), breaks = 1:10) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "eigenvalues", y = "variance explained (%)")

ggsave(paste0(folder_data, "temp/18-01-snp_PCs.png"), p, width = 4, height = 3)

# 2. PCA based on SNPs ----
set.seed(1)
p1 <- sc %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = rhizobia_population), shape = 21, size = 3, stroke = 1, position = position_jitter(width = 10, height = 10)) +
    scale_color_manual(values = rhizobia_population_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 variance = ", pcs[1], "%"), y = paste0("PC2 variance = ", pcs[2], "%"))

}

