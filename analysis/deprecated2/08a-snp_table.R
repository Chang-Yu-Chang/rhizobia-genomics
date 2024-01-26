
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

snps <- read_delim("~/Dropbox/lab/local-adaptation/data/genomics/popgen/read_em1021/snippy/core.tab", delim = "\t")
snps <- clean_names(snps)
dim(snps)
table(snps$chr)

snps <- read_delim("~/Dropbox/lab/local-adaptation/data/genomics/popgen/read_wsm419/snippy/core.tab", delim = "\t")
snps <- clean_names(snps)
dim(snps)
table(snps$chr)


# isolates_mash <- read_csv(paste0(folder_data, "temp/14-isolates_mash.csv")) 

# isolates_mash %>%
#     select(genome_id, species_name) %>%
#     arrange(species_name) %>%
#     view
    
