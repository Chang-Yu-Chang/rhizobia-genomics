#' This script generates Table S1

renv::load()
library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
en <- read_csv(paste0(folder_data, "raw/ensifer_ncbi.csv"), col_names = F)
en <- en[,1:3]
colnames(en) <- c("Accession", "Species", "Strain")
en <- en %>%
    mutate(Species = str_replace(Species, "E.", "Ensifer")) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

# Flextable
ft <- flextable(en) %>%
    autofit() %>%
    merge_v(j = "Species") %>%
    valign(valign = "top") %>%
    style(j = "Species", pr_t = fp_text_default(italic = T))

save_as_image(ft, path = here::here("plots/TabS1.png"))
