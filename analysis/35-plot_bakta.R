#' This script plot the annotation by bakta

library(tidyverse)
library(cowplot)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# Read the contig information
g_contigs <- read_csv(paste0(folder_data, "temp/34-g_contigs.csv"), show_col_types = F)
ann_bak <- read_tsv(paste0(folder_data, "temp/plasmidsaurus/Chang_Q5C_2/05-bakta/consensus.tsv"), skip = 2, show_col_types = F)
ann_bak <- ann_bak %>% clean_names()
#
g_contigs
ann_bak %>%
    group_by(number_sequence_id) %>%
    summarize(n())


ann_bak %>%
    #filter(str_detect(product, "plasmid"))
    filter(str_detect(gene, "fix")) %>%
    view


"HALT ANAYLSIS ON bakta because the g1 contig2 is the largeest in medaka but in bakta it's contig3. Have to rerun bakta"
