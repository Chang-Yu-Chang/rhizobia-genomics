#' Clean the list of symbiosis genes

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

tb <- read_csv("~/Downloads/symbiosis_genes.csv", col_names = F)

# Long format
tb <- tb %>%
    mutate(temp = 1:n()) %>%
    pivot_longer(cols = -temp) %>%
    drop_na(value) %>%
    select(gene = value)

nrow(tb) # 561 genes

# Clean names
tb <- tb %>%
    mutate(
        gene_id = 1:n(),
        gene = str_remove(gene, "\\[\\d+\\]"), # remove reference
        temp = str_extract(tb$gene, "\\(([^)]+)\\)") %>%str_remove_all("\\(|\\)"), # extract parentheses
        gene = str_remove(gene, "\\(([^)]+)\\)")
    ) %>%
    pivot_longer(-gene_id) %>%
    select(-name, gene = value) %>%
    drop_na(gene) %>%
    distinct(gene)


write_csv(tb, paste0(folder_data, "genomics_analysis/gene_content/symbiosis_genes.csv"))
