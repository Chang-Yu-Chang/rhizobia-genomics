#' This script analyse the raw growth curve data

library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

rhizobia_pilot <- readxl::read_xlsx(paste0(folder_data, "raw/growth_curve/rhizobia_pilot.xlsx"))
rhizobia_pilot_plate <- read_csv(paste0(folder_data, "raw/growth_curve/rhizobia_pilot_plate.csv"))
clean_well_names <- function (x) {
    y <- paste0(str_sub(x, 1, 1), str_sub(x, 2, 3) %>% as.numeric %>% sprintf("%02d", .))
    return(y)
}

# Tidy up
rhizobia_pilot <- rhizobia_pilot %>%
    # Calculate the time interval in minutes
    mutate(Time = (rhizobia_pilot$Time[1] %--% Time)/dminutes(1)) %>%
    #mutate(Time = hour(rhizobia_pilot$Time) + minute(rhizobia_pilot$Time)/60) %>%
    select(-`T° 600`) %>%
    pivot_longer(cols = -Time, names_to = "Well", values_to = "OD600") %>%
    mutate(Well = clean_well_names(Well)) %>%
    left_join(rhizobia_pilot_plate)


# 1. Raw data
p <- rhizobia_pilot %>%
    ggplot() +
    geom_line(aes(x = Time, y = OD600, color = Strain, group = Well)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04-01-rhizobia_pilot_raw.png"), p, width = 6, height = 5)

# 2. Curve by well
p <- rhizobia_pilot %>%
    mutate(row = str_sub(Well, 1, 1), col = str_sub(Well, 2, 3)) %>%
    ggplot() +
    geom_line(aes(x = Time, y = OD600, color = Strain), linewidth = 1) +
    theme_light() +
    facet_grid(row ~ col) +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04-02-rhizobia_pilot_raw_grid.png"), p, width = 20, height = 15)

# 3. Average OD by strain
rhizobia_pilot %>%
    group_by(Time, Strain) %>%
    summarize(mean_OD600 = mean(OD600), sd_OD600 = sd(OD600)) %>%
    ggplot() +
    geom_line(aes(x = Time, y = mean_OD600, color = Strain)) +
    geom_ribbon(aes(x = Time, ymin = mean_OD600 - sd_OD600, ymax = mean_OD600 + sd_OD600, fill = Strain), alpha = 0.2) +
    theme_classic() +
    theme() +
    guides() +
    labs()



# rhizobia_pilot_plate <- tibble(Well = paste0(rep(LETTERS[1:8], 12), rep(sprintf("%02d", 1:12), each = 8)),
#        Strain = c(rep(c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2"), each = 8), rep("blank", 6*8)))
#
# write_csv(rhizobia_pilot_plate, paste0(folder_data, "raw/growth_curve/rhizobia_pilot_plate.csv"))





















