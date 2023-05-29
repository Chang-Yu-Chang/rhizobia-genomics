#' This script is used for analyzing the pilot test on growth curve
#' The experiment was done on 20230517
#' Read the README.txt in raw/growth_curve/ for details

library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

pilot <- readxl::read_xlsx(paste0(folder_data, "raw/growth_curve/pilot.xlsx"))
pilot_plate <- read_csv(paste0(folder_data, "raw/growth_curve/pilot_plate.csv"))
clean_well_names <- function (x) {
    y <- paste0(str_sub(x, 1, 1), str_sub(x, 2, 3) %>% as.numeric %>% sprintf("%02d", .))
    return(y)
}

# Tidy up
pilot <- pilot %>%
    # Calculate the time interval in minutes
    mutate(Time = (pilot$Time[1] %--% Time)/dminutes(1)) %>%
    #mutate(Time = hour(pilot$Time) + minute(pilot$Time)/60) %>%
    select(-`T° 600`) %>%
    pivot_longer(cols = -Time, names_to = "Well", values_to = "OD600") %>%
    mutate(Well = clean_well_names(Well)) %>%
    left_join(pilot_plate)


# 1. Raw data
p <- pilot %>%
    ggplot() +
    geom_line(aes(x = Time, y = OD600, color = Strain, group = Well)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04-01-gc_pilot_raw.png"), p, width = 6, height = 5)

# 2. Curve by well
p <- pilot %>%
    mutate(row = str_sub(Well, 1, 1), col = str_sub(Well, 2, 3)) %>%
    ggplot() +
    geom_line(aes(x = Time, y = OD600, color = Strain), linewidth = 1) +
    theme_light() +
    facet_grid(row ~ col) +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04-02-gc_pilot_raw_grid.png"), p, width = 20, height = 15)

# 3. Average OD by strain
p <- pilot %>%
    group_by(Time, Strain) %>%
    summarize(mean_OD600 = mean(OD600), sd_OD600 = sd(OD600)) %>%
    ggplot() +
    geom_line(aes(x = Time, y = mean_OD600, color = Strain)) +
    geom_ribbon(aes(x = Time, ymin = mean_OD600 - sd_OD600, ymax = mean_OD600 + sd_OD600, fill = Strain), alpha = 0.2) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04-03-gc_pilot_mean.png"), p, width = 6, height = 5)


















