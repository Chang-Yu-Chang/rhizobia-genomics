#'

library(tidyverse)
library(janitor)
# library(cowplot)
# library(ggh4x)
source(here::here("metadata.R"))

tb1 <- read_csv(paste0(folder_data, "variance/March_DZA_final.csv"))
tb2 <- read_csv(paste0(folder_data, "variance/Sept_DZA_final.csv"))
tb3 <- read_csv(paste0(folder_data, "variance/May_A17_final.csv"))
tb4 <- read_csv(paste0(folder_data, "variance/Nov_A17_final.csv"))
strains <- read_csv(paste0(folder_data, "variance/strain_metadata.csv")) %>% clean_names()

tb <- bind_rows(tb1, tb2, tb3, tb4) %>%
    clean_names() %>%
    arrange(exp, plant_id, strain_id)

sc <- function (x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))

tb_strains <- tb %>%
    group_by(strain_id, line) %>%
    filter(strain_id != "control_ex") %>%
    mutate(across(c(nod, nod_weight, shoot, root), sc)) %>%
    summarize(across(c(nod, nod_weight, shoot, root), function(x) mean(x, na.rm = T)), n())

tb_strains %>%
    group_by(line) %>%
    summarize(across(c(nod, nod_weight, shoot, root), function (x) var(x, na.rm = T)))

tb_strains %>%
    ggplot() +
    geom_histogram(aes(x = root, fill = line), position = "identity", alpha = .5) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()
