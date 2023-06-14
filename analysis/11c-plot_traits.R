#' This script plots the root traits

library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

# growth curves
gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)

# common garden
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/11-treatments_long.csv"), show_col_types = F)
treatments_scaled <- read_csv(paste0(folder_data, "temp/11-treatments_scaled.csv"), show_col_types = F)
treatments_scaled_long <- read_csv(paste0(folder_data, "temp/11-treatments_scaled_long.csv"), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(rhizobia = ExpID)

# 1.
p <- gc_summ %>%
    mutate(strain = factor(strain, list_strains)) %>%
    filter(strain == "H1M1R1") %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs), color = "maroon") +
    geom_ribbon(aes(x = t, ymin = mean_abs - sd_abs, ymax = mean_abs + sd_abs), alpha = 0.2, fill = "maroon") +
    scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, 12), expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.04, 0.6), expand = c(0,0)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "time (hrs)", y = expression(OD[600]))

ggsave(paste0(folder_data, "temp/04a-05-gc_mean_example.png"), p, width = 4, height = 3)
