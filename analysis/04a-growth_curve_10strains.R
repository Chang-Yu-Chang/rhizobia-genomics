#' This script is used for analyzing the growth curve
#' The experiment was done on 20230526
#' Read the README.txt in raw/growth_curve/ for details

library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

gc <- read_csv(paste0(folder_data, "raw/growth_curve/rhizobia_growth_curve.csv"))
gc_plate <- read_csv(paste0(folder_data, "raw/growth_curve/gc_plate.csv"))
clean_well_names <- function (x) {
    y <- paste0(str_sub(x, 1, 1), str_sub(x, 2, 3) %>% as.numeric %>% sprintf("%02d", .))
    return(y)
}

# Tidy up
gc <- gc %>%
    # Calculate the time interval in minutes
    mutate(time_step = as.numeric(difftime(lead(Time), Time[1], units = "hours"))) %>%
    #mutate(Time = hour(gc$Time) + minute(gc$Time)/60) %>%
    select(-Time, -`T° 600`) %>%
    pivot_longer(cols = -time_step, names_to = "well", values_to = "od600") %>%
    mutate(well = clean_well_names(well)) %>%
    left_join(gc_plate)


# 1. Raw data
p <- gc %>%
    ggplot() +
    geom_line(aes(x = time_step, y = od600, color = strain, group = well)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-01-gc_raw.png"), p, width = 6, height = 5)

# 2. Curve by well
p <- gc %>%
    mutate(row = str_sub(well, 1, 1), col = str_sub(well, 2, 3)) %>%
    ggplot() +
    geom_line(aes(x = time_step, y = od600, color = strain), linewidth = 1) +
    theme_light() +
    facet_grid(row ~ col) +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-02-gc_raw_grid.png"), p, width = 20, height = 15)

# 3. Average OD by strain
gc_summ <- gc %>%
    group_by(time_step, strain) %>%
    summarize(mean_od600 = mean(od600), sd_od600 = sd(od600))
p <- gc_summ %>%
    filter(strain %in% c("H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1", "L4M2R2", "L4M3R3", "L4M4R1", "blank")) %>%
    mutate(strain = factor(strain, c("H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1", "L4M2R2", "L4M3R3", "L4M4R1", "blank"))) %>%
    ggplot() +
    geom_line(aes(x = time_step, y = mean_od600, color = strain)) +
    geom_ribbon(aes(x = time_step, ymin = mean_od600 - sd_od600, ymax = mean_od600 + sd_od600, fill = strain), alpha = 0.2) +
    facet_wrap(~strain, ncol = 5) +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-03-gc_mean.png"), p, width = 10, height = 6)

# 4. plot the growth curve fit
gc.prm.stat <- read_csv(file=paste0(folder_data, 'temp/04c-gc_prm_summ.csv'))
gc.prm <- read_csv(file=paste0(folder_data, 'temp/04c-gc_prm.csv'))

gc.prm.stat %>%
    filter(strain %in% c("H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1", "L4M2R2", "L4M3R3", "L4M4R1", "blank")) %>%
    ggplot() +
    geom_point(aes(x = strain, y = r)) +
    theme_classic() +
    theme() +
    guides() +
    labs()


p <- gc.prm %>%
    filter(strain %in% c("H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1", "L4M2R2", "L4M3R3", "L4M4R1", "blank")) %>%
    mutate(strain = factor(strain, rev(c("H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1", "L4M2R2", "L4M3R3", "L4M4R1", "blank")))) %>%
    pivot_longer(cols = c(r, t.r, lag, maxOD.fit), names_to = "gc_trait") %>%
    ggplot() +
    geom_boxplot(aes(x = strain, y = value), outlier.size = 0) +
    geom_point(aes(x = strain, y = value), shape = 21) +
    facet_wrap(gc_trait ~., scales = "free", ncol = 1, strip.position = "right") +
    coord_flip() +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/04a-04-gc_trait.png"), p, width = 6, height = 10)











