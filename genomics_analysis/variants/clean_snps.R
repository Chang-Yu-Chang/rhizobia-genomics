#' This script cleans the SNPs data

renv::load()
library(tidyverse)
library(proxy) # For computing jaccard distance
library(vcfR) # for handling VCF
library(poppr) # for pop gen analysis
library(vegan) # for analysis of variance using distance matrices
source(here::here("metadata.R"))

isolates_contigs <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_contigs.csv"))
snp <- read_table(paste0(folder_data, "genomics/variants/em1021/core.tab"))
gs <- read_table(paste0(folder_data, "genomics/variants/em1021/core.txt"))
vcf <- read.vcfR(paste0(folder_data, "genomics/variants/em1021/core.vcf"))

#
vcfR::getCHROM(vcf) %>% table() # SNPs on chromosome
gl <- vcfR2genlight(vcf) # VCF to genlight
snps <- tab(gl); dim(snps) # Number of SNPs 18450

# 1. SNPs
## Test VA populations: high vs low elevation
isolates1 <- isolates %>% filter(population == "VA") %>% left_join(isolates_contigs)
m <- snps[match(isolates1$genome_id, gl@ind.names),]; dim(m) # Filter the snps table for focal genomes
#m <- m[,apply(m, 2, sum) != 0]; dim(m) # 18338 snps
dist_m1 <- dist(m, "jaccard")
permanova1 <- adonis2(dist_m1 ~ site_group, data = isolates1)
permanova1 # no
pcoa1 <- cmdscale(dist_m1, eig = T)
eigs1 <- round(pcoa1$eig / sum(pcoa1$eig)*100, 2)

## Test PA populations: urban vs suburban elevation
isolates2 <- isolates %>% filter(population == "PA") %>% left_join(isolates_contigs)
m <- snps[match(isolates2$genome_id, gl@ind.names),]; dim(m) # Filter the snps table for focal genomes
m <- m[,apply(m, 2, sum) != 0]; dim(m) # 323 snps
dist_m2 <- dist(m, "jaccard")
permanova2 <- adonis2(dist_m2 ~ site_group, data = isolates2)
permanova2 # no
pcoa2 <- cmdscale(dist_m2, eig = T)
eigs2 <- round(pcoa2$eig / sum(pcoa2$eig)*100, 2)


save(isolates1, isolates2, pcoa1, pcoa2, eigs1, eigs2, file = paste0(folder_data, "genomics_analysis/variants/snps.rdata"))
