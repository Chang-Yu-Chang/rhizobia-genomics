#' This script plots the MK test results

library(tidyverse)
library(janitor)
library(ggsci)
library(cowplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))


set_name = "elev_med"
ref = "ngr234"
mktests1 <- read_csv(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/mktests.csv"))
tt1 <- read_gpas(set_name)

set_name = "urbn_mel"
ref = "ngr234"
mktests2 <- read_csv(paste0(folder_data, "genomics_analysis/mktest/", set_name, "/", ref, "/mktests.csv"))
tt2 <- read_gpas(set_name)


mm1 <- mktests1 %>%
    left_join(select(tt1$cleaned_gene_names, gene, from))
mm2 <- mktests2 %>%
    left_join(select(tt2$cleaned_gene_names, gene, from))
