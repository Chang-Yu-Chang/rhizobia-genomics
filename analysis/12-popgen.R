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

core_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy/core/core.vcf"))
medicae_vcf <- read.vcfR(paste0(folder_data, "temp/anvio/06-alignment/snippy_medicae/core/core.vcf"))

vcfR::getCHROM(medicae_vcf) %>% table() %>% sum()
