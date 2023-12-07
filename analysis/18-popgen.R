#' This script runs population genetics analysis

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(ggsci)
    library(vcfR) # for handling VCF
    library(poppr) # for pop gen analysis
    source(here::here("analysis/00-metadata.R"))
})

snp_usda <- read.vcfR(paste0(folder_data, "genomics/popgen/snippy_usda1106/snippy_usda1106.vcf"))
snp_wsm <- read.vcfR(paste0(folder_data, "genomics/popgen/snippy_wsm419/snippy_wsm419.vcf"))

# 1. descriptions ----
# number of snps on each contig
vcfR::getCHROM(snp_usda) %>% table()
vcfR::getCHROM(snp_wsm) %>% table()


# Convert the vcf to a genclone object
geni_usda <- vcfR2genind(snp_usda)
# core_genind <- as.genclone(core_genind)
isolates_ensifer <- read_csv(paste0(folder_data, "temp/02-isolates_rhizo.csv"), show_col_types = F) %>%
    filter(genus == "Ensifer") %>%
    mutate(genome_name = str_replace(genome_id, "g", "Chang_Q5C_")) %>%
    select(genome_name, everything()) %>%
    mutate(genome_name = factor(genome_name, rownames(core_genind$tab))) %>%
    arrange(genome_name) %>%
    drop_na()
other(core_genind)$taxa <- isolates_ensifer
