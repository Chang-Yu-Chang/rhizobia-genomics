#' This script generates the table of ncbi accessions

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Table
en <- read_csv(paste0(folder_data, "raw/ensifer_ncbi.csv"), col_names = F)
en <- en[,1:3]
colnames(en) <- c("Accession", "Species", "Strain")
en <- en %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

# Flextable
ft <- flextable(en) %>%
    autofit() %>%
    merge_v(j = "Species") %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = c(1:4,8:11,13,17, 19), j = 1:4) %>%
    align(part = "all", align = "center") %>%
    valign(j = "Species", valign = "top") %>%
    style(j = "Species", pr_t = fp_text_default(italic = T)) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS2.html"))
#save_as_image(ft, path = here::here("plots/TabS2.png"))
