#' This script uses DAYMET https://daymet.ornl.gov/ database to extract the
#' climate for our sampling sites given the coordinate

renv::load()
library(daymetr)
library(tidyverse)
library(janitor)
library(RColorBrewer)
library(elevatr) # for getting elevation
library(sf)
source(here::here("analysis/00-metadata.R"))


dms2dec <- function(x) {
    #' convert degree-minute-second coordinate to decimal
    str_split(x, " ") %>%
        sapply(function(xx) {
            z <- as.numeric(xx[1:3])
            s <- ifelse(xx[4] %in% c("N", "E"), 1, -1)
            s*(z[1] + z[2]/60 + z[3]/3600)
        }) %>%
        return
}
sites_mlbs <- readxl::read_excel(paste0(folder_data, "raw/High and low elevation sites_Sept 2022.xlsx"))
sites_phila <- read_csv(paste0(folder_data, "raw/sites_phila.csv"), show_col_types = F) %>%
    drop_na() %>%
    mutate(description = tolower(description) %>% str_remove("'")) %>%
    separate(col = coord, into = c("latitude_dec", "longitude_dec"), sep = ",\\s", convert = T)



# 0. Clean up the site coordinates ----
# mlbs
sites_mlbs <- sites_mlbs %>%
    clean_names() %>%
    mutate(latitude_dec = dms2dec(latitude_degrees_minutes_seconds)) %>%
    mutate(longitude_dec = dms2dec(longitude_degrees_minutes_seconds)) %>%
    mutate(elevation_m = elevation_feet / 3.28084) %>%
    mutate(site_group = str_sub(site,1 , 1)) %>%
    select(site_group, site, latitude_dec, longitude_dec, elevation_m)

# phila
temp <- sites_phila %>%
    st_as_sf(coords = c("longitude_dec", "latitude_dec"), crs = 4326) %>%
    get_elev_point()

sites_phila <- sites_phila %>%
    mutate(elevation_m = temp$elevation) %>%
    select(site_group, site, latitude_dec, longitude_dec, elevation_m)

# bind rows
sites <- bind_rows(mutate(sites_mlbs, population = "MLBS"), mutate(sites_phila, population = "Phila")) %>%
    select(population, everything())

write_csv(sites, paste0(folder_data, "temp/22-sites.csv"))

# 1. Extract the climatology data from Daymet ----
list_dm <- rep(list(NA), nrow(sites))
for (i in 1:nrow(sites)) {
    dm <- download_daymet(site = sites$site[i],
                          lat = sites$latitude_dec[i],
                          lon = sites$longitude_dec[i],
                          start = 2022,
                          end = 2022,
                          internal = TRUE)
    list_dm[[i]] <- as_tibble(dm$data) %>% mutate(site = sites$site[i], population = sites$population[i])

}

dml <- bind_rows(list_dm) %>%
    mutate(site_group = str_sub(site,1 , 1)) %>%
    clean_names() %>%
    mutate(ydate = strptime(paste("2022", yday), format="%Y %j")) %>%
    mutate(ymonth = month(ydate))


write_csv(dml, paste0(folder_data, "temp/22-dml.csv"))



# 2. Resample the difference between a H and a L temperature ---
## tmax
set.seed(1)
n_resample <- 100
list_diff_tmax <- rep(list(tibble(resample = 1:n_resample)), 365)
for (d in 1:365) {
    tmax_H <- dml %>% filter(site_group == "H", yday == d) %>% pull(tmax_deg_c)
    tmax_L <- dml %>% filter(site_group == "L", yday == d) %>% pull(tmax_deg_c)
    sampled_tmax_H <- sample(tmax_H, n_resample, replace = T)
    sampled_tmax_L <- sample(tmax_L, n_resample, replace = T)
    list_diff_tmax[[d]]$yday <- d
    list_diff_tmax[[d]]$diff_tmax <- sampled_tmax_L - sampled_tmax_H
    cat("\t", d)
}

diff_tmax <- bind_rows(list_diff_tmax)
write_csv(diff_tmax, paste0(folder_data, "temp/22-diff_tmax.csv"))

## tmin
set.seed(1)
n_resample <- 100
list_diff_tmin <- rep(list(tibble(resample = 1:n_resample)), 365)
for (d in 1:365) {
    tmin_L <- dml %>% filter(site_group == "L", yday == d) %>% pull(tmin_deg_c)
    tmin_H <- dml %>% filter(site_group == "H", yday == d) %>% pull(tmin_deg_c)
    sampled_tmin_H <- sample(tmin_H, n_resample, replace = T)
    sampled_tmin_L <- sample(tmin_L, n_resample, replace = T)
    list_diff_tmin[[d]]$yday <- d
    list_diff_tmin[[d]]$diff_tmin <- sampled_tmin_L - sampled_tmin_H
    cat("\t", d)
}

diff_tmin <- bind_rows(list_diff_tmin)
write_csv(diff_tmin, paste0(folder_data, "temp/22-diff_tmin.csv"))



# 3. Shade for month ----
tb_season <- dml %>%
    distinct(yday, ydate, ymonth) %>%
    mutate(season = case_when(
        ymonth %in% 1:3 ~ "spring",
        ymonth %in% 4:6 ~ "summer",
        ymonth %in% 7:9 ~ "fall",
        ymonth %in% 10:12 ~ "winter",
    )) %>%
    group_by(season) %>%
    filter(ydate == max(ydate) | ydate == min(ydate)) %>%
    mutate(temp = c("start", "end")) %>%
    select(-ydate, -ymonth) %>%
    pivot_wider(names_from = temp, values_from = yday) %>%
    ungroup()

tb_summer <- dml %>%
    distinct(yday, ydate, ymonth) %>%
    mutate(season = case_when(
        ymonth %in% 5:10 ~ "summer",
        T ~ NA
    )) %>%
    group_by(season) %>%
    filter(ydate == max(ydate) | ydate == min(ydate)) %>%
    mutate(temp = c("start", "end")) %>%
    drop_na(season) %>%
    select(-ydate, -ymonth) %>%
    pivot_wider(names_from = temp, values_from = yday) %>%
    ungroup()

tb_month <- dml %>%
    distinct(yday, ydate, ymonth) %>%
    group_by(ymonth) %>%
    filter(ydate == max(ydate) | ydate == min(ydate)) %>%
    mutate(temp = c("start", "end")) %>%
    select(-ydate) %>%
    pivot_wider(names_from = temp, values_from = yday) %>%
    mutate(ymonth = factor(ymonth)) %>%
    ungroup()

write_csv(tb_summer, paste0(folder_data, "temp/22-tb_summer.csv"))
write_csv(tb_season, paste0(folder_data, "temp/22-tb_season.csv"))
write_csv(tb_month, paste0(folder_data, "temp/22-tb_month.csv"))
