#' This script plots the MK test results

library(tidyverse)
library(janitor)
library(ggsci)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

read_mktests <- function (set_name, ref) {
    # set_name = "elev_med"
    # ref = "ngr234"
    mktests <- read_csv(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/mktests.csv"))
    blast <- read_table(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/blast_results.txt"), col_names = F) %>%
        select(1:4) %>%
        setNames(c("gene", "refseq", "pident", "length"))
    return(list(mktests = mktests, blast = blast))
}

tt <- read_gpas("elev_med")
mm <- read_mktests("elev_med", "ngr234")
# tt2 <- read_gpas("urbn_mel")
# mm2 <- read_mktests("urbn_mel", "ngr234")

tt$list_sccg %>%
    filter(str_detect(gene, "aac"))

mm1$mktests %>%
    left_join(mm1$blast)

mm1$blast %>%
    distinct(gene)

mm1$blast %>%
    group_by(gene) %>%
    count() %>%
    filter(n!=1)
    distinct(gene)
    view


mm1 <- mktests1 %>%
    left_join(select(tt1$cleaned_gene_names, gene, from))
mm2 <- mktests2 %>%
    left_join(select(tt2$cleaned_gene_names, gene, from))

max(mm1$mktest_alpha)

mm1 %>%
    drop_na(from) %>%
    filter(from %in% mm2$from) %>%
    left_join(select(mm2, from, gene2 = gene, mktest_alpha2 = mktest_alpha), relationship = "many-to-many") %>%
    ggplot() +
    geom_point(aes(x = mktest_alpha, y = mktest_alpha2)) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_bw() +
    theme() +
    guides() +
    labs()

