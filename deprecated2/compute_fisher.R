#' This script computes Fisher statistic for each gene presence/absence

renv::load()
library(tidyverse)
library(broom)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
genomes <- read_csv(paste0(folder_data, "mapping/genomes.csv"))
isolates_tax <- read_csv(paste0(folder_data, "genomics_analysis/taxonomy/isolates_tax.csv"))
gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gpa.csv"))
gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/gene_order.csv"))
#which(gene_order == "fcl_2")

# For all genes in the pangenome
tb_fis <- tibble(gene = gpa$gene, fis = NA)

for (i in 1:nrow(gpa)) {
    if (i %% 100 == 0) print(i)

    tball <- gpa[i,] %>%
        pivot_longer(-gene, names_to = "genome_id", values_to = "presence") %>%
        left_join(isolates, by = join_by(genome_id)) %>%
        mutate(presence = factor(presence, c(0, 1)))

    tball_elev <- tball %>% filter(population == "VA")
    tball_urba <- tball %>% filter(population == "PA")

    fis1 <- table(tball_elev$site_group, tball_elev$presence) %>%
        # Haldane-Anscombe correction
        `+`(1) %>%
        fisher.test() %>%
        tidy()

    fis2 <- table(tball_urba$site_group, tball_urba$presence) %>%
        `+`(1) %>%
        fisher.test() %>%
        tidy()


    tb_fis$fis[i] <- list(bind_rows(
        tibble(gradient = "elevation") %>% bind_cols(fis1),
        tibble(gradient = "urbanization") %>% bind_cols(fis2)
    ))
}

gpa_fis <- tb_fis %>% unnest(fis)

write_csv(gpa_fis, paste0(folder_data, "genomics_analysis/fst/gpa_fis.csv"))


