#' This script runs population genetics analysis

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(ggsci)
    library(vcfR) # for handling VCF
    #library(VariantAnnotation)
    library(poppr) # for pop gen analysis
    source(here::here("analysis/00-metadata.R"))
})

isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F)
isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"), show_col_types = F)
snp_usda <- read_table(paste0(folder_data, "genomics/popgen/snippy_usda1106/snippy_usda1106.tab"), show_col_types = F)
vcf_usda <- read.vcfR(paste0(folder_data, "genomics/popgen/snippy_usda1106/snippy_usda1106.vcf"))


# VCF to genlight
gl_usda <- vcfR2genlight(vcf_usda)

#pop(gl_usda) <- isolates_mapping$rhizobia_site
#ploidy(gl_usda)
#gl_usda@ind.names <- isolates_mapping$genome_id


# PCA
pca_usda <- glPca(gl_usda, nf = 2)

# convert scores of vcf.pca into a tibble
sc_usda <- as_tibble(pca_usda$scores)




# add the country data into a column of vcf.pca.scores tibble
vcf.pca.scores$country <- metadata$country


# We will also determine the variance each PC contributes the data, which will help us understand potential drivers of patterns in our dataset. Lets plot the eigenvectors to try an understand this a bit more.

barplot(100 * vcf.pca$eig / sum(vcf.pca$eig), col="green")
title(ylab = "Percent of variance explained")
title(xlab = "Eigenvalues")



vcfR::getCHROM(vcf_usda) %>% table()

gi_usda <- vcfR2genind(vcf_usda)
class(gi_usda)
gc_usda <- as.genclone(gi_usda)
class(gc_usda)

gc_usda$pop <- factor(isolates_mapping$rhizobia_site)


glPca(gc_usda)



popsub(gc_usda, sublist = "urban", drop = T)

amova(gc_usda)






if (FALSE) {
    vcf_usda <- readVcf(paste0(folder_data, "genomics/popgen/snippy_usda1106/snippy_usda1106.vcf"))
vcf_wsm <- read.vcfR(paste0(folder_data, "genomics/popgen/snippy_wsm419/snippy_wsm419.vcf"))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"), show_col_types = F)
isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"), show_col_types = F)

# 1. VCF description ----
# number of snps on each contig
vcfR::getCHROM(vcf_usda) %>% table()
vcfR::getCHROM(vcf_wsm) %>% table()

# extract only the chromosome

as_tibble(vcf_usda@fix) %>%
    distinct(POS)
snp_usda <- as_tibble(vcf_usda@gt, .name_repair = "minimal") %>%
    setNames(c("snp", paste0("g", 1:32)))



snp_usda %>%
    clean_names() %>%
    slice(1:10000) %>%
    pivot_longer(cols = -snp) %>%
    filter(value != 1)


# 1. descriptions ----

# Convert the vcf to a genclone object
gi_usda <- vcfR2genind(snp_usda)
class(gi_usda)
gc_usda <- as.genclone(gi_usda)
class(gc_usda)

gc_usda$pop <- factor(isolates_mapping$rhizobia_site)


poppr(gc_usda)


# mll(geni_usda)
#
# mll(monpop) <- "original"
# diversity_stats(mll(monpop))
#mlg.filter(x, distance = xdis) <- 1 + .Machine$double.eps^0.5



# Remove uninformative loci
#geni_usda_cut <- informloci(geni_usda, cutoff = 2/nInd(geni_usda), MAF = 0.01, quiet = FALSE)
popNames(core_genind)


genlight_usda <- dartR::gi2gl(geni_usda)












pca_usda <- glPca(genlight_usda, nf = 2, center = T, scale = T)
eig_usda <- sort(pca_usda$eig/sum(pca_usda$eig), decreasing = T)

# Merge the tb
isolates_pca <- pca_usda$scores %>%
    as_tibble() %>%
    mutate(PC1 = PC1/100, PC2 = PC2/100) %>%
    mutate(genome_name = genlight_usda$ind.names) %>%
    left_join(isolates_ensifer)

set.seed(1)
p1 <- isolates_pca %>%
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = site), shape = 21, size = 3, stroke = 1.5, position = position_jitter(width = 0.02, height = 0.02)) +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        #legend.background = element_rect(fill = "white", color = "black"),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.2, linetype = 2),
        panel.border = element_rect(fill = NA, color = "black")
    ) +
    guides() +
    labs(x = paste0("PC1 (", round(eig[1]*100,1), "%)"),
         y = paste0("PC2 (", round(eig[2]*100,1), "%)"))
}
