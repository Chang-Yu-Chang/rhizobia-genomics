#'

library(tidyverse)
library(cowplot)
source(here::here("metadata.R"))
library(terra)

fi <- rast("~/Downloads/NLCD_rjyHFk4exYhJdCwbUUXp/Annual_NLCD_FctImp_2022_CU_C1V0_rjyHFk4exYhJdCwbUUXp.tiff")
plot(fi)

fi_df <- project(fi, "EPSG:4326") %>%
    terra::as.data.frame(xy = TRUE, na.rm = TRUE) %>%
    as_tibble
names(fi_df)[3] <- "fraimp"

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv")) %>% filter(population != "mid elevation") %>%
    # Use sites where the isolates were from
    filter(site %in% isolates$site)
sites_center <- sites %>%
    group_by(gradient) %>%
    summarize(lat_mean = mean(latitude_dec), lon_mean = mean(longitude_dec), lat_max = max(latitude_dec), lat_min = min(latitude_dec), lon_max = max(longitude_dec), lon_min = min(longitude_dec)) %>%
    mutate(hjust = c(0.8, 0.8)) %>%
    mutate(width_max = max(c((lon_max-lon_min)*1.8, (lat_max-lat_min)*1.8)))


p <- fi_df %>%
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = fraimp)) +
    scale_fill_gradient2(low = "#a642bf", high = "#0cc45f", mid = "snow", midpoint = 50, name = "frac. impervious surface", breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
    geom_point(data = filter(sites, gradient == "urbanization"), aes(x = longitude_dec, y = latitude_dec, color = population), fill = "black", size = 2, shape = 21, stroke = 2) +
    scale_x_continuous(limits = c(sites_center$lon_mean[2] - sites_center$width_max[2]/2, sites_center$lon_mean[2] + sites_center$width_max[2]/2), expand = c(0,0)) +
    scale_y_continuous(limits = c(sites_center$lat_mean[2] - sites_center$width_max[2]/2, sites_center$lat_mean[2] + sites_center$width_max[2]/2), expand = c(0,0)) +
    coord_fixed() +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_phenotypes, "sites/fracimp.png"), p, width = 9, height = 7)
