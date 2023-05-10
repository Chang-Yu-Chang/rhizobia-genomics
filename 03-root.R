#' This script analyses the hand measured phenotypes: dry weight, nodule count, root weight, nodule weight

library(tidyverse)
library(broom)
library(janitor)
source(here::here("analysis/00-metadata.R"))

#
treatments <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(dry_weight = `DryWeight (mg)`) %>%
    clean_names()

features <- read_csv(paste0(folder_data, "raw/rhizobia/05-root_architecture/features.csv"), show_col_types = F) %>%
    clean_names() %>%
    mutate(id = str_replace(file_name, ".png", "") %>% as.numeric())


treatments <- treatments %>% left_join(features)

#write_csv(treatments, paste0(folder_data, "temp/03-treatments.csv"))

#
treatments <- read_csv(paste0(folder_data, "temp/03-treatments.csv"), show_col_types = F)
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), unique(treatments$rhizobia))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")


#
treatments %>%
    ggplot() +
    geom_point(aes(x = dry_weight, y = network_area_px2)) +
    theme_classic() +
    theme() +
    guides() +
    labs()
