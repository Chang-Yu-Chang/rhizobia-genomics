#' This script plots the MK test results

library(tidyverse)
library(janitor)
library(ggsci)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

read_mktests <- function (set_name, ref) {
    #' Read the blast and MK test results
    mktests <- read_csv(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/mktests.csv"))
    blast <- read_table(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/blast_results.txt"), col_names = F) %>%
        select(1:4) %>%
        setNames(c("gene", "refseq", "pident", "length"))
    return(list(mktests = mktests, blast = blast))
}

# Join the two gradients
tt1 <- read_gpas("elev_med")
mm1 <- read_mktests("elev_med", "ngr234")
tt2 <- read_gpas("urbn_mel")
mm2 <- read_mktests("urbn_mel", "ngr234")

mk1 <- mm1$mktests %>%
    left_join(select(tt1$cleaned_gene_names, gene, from)) %>%
    drop_na(from)

mk2 <- mm2$mktests %>%
    left_join(select(tt2$cleaned_gene_names, gene, from)) %>%
    drop_na(from)


mks <- mk1 %>%
    left_join(select(mk2, from, gene2 = gene, mktest_alpha2 = mktest_alpha), relationship = "many-to-many")

cor.test(mks$mktest_alpha, mks$mktest_alpha2)
#lm(mktest_alpha2 ~ mktest_alpha, data = mks) %>%  summary()

mks <-
    ggplot() +
    geom_point(aes(x = mktest_alpha, y = mktest_alpha2)) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_bw() +
    theme() +
    guides() +
    labs()

