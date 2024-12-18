#' This script uses DAYMET https://daymet.ornl.gov/ database to extract the climate data for our sampling sites given the coordinates
#' 0. Clean the site coordinate format and obtain elevation data
#' 1. Extract the year 2022 climatology data from Daymet
#' 2. Bootstrap the temp difference between populations
#' 3. Get the day for plotting the month shades

library(tidyverse)
library(janitor)
library(daymetr) # for extracting climate data
library(elevatr) # for getting elevation
library(sf) # for handing simple features
Sys.setenv(PROJ_LIB = "/opt/homebrew/Cellar/proj/9.5.0/share/proj") # for crs
source(here::here("metadata.R"))

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
sites_mlbs <- readxl::read_excel(paste0(folder_data, "raw/plants/High and low elevation sites_Sept 2022.xlsx"))
sites_phila <- read_csv(paste0(folder_data, "raw/plants/sites_phila.csv"), show_col_types = F) %>%
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
    mutate(site = str_remove(site, "-")) %>%
    mutate(population = case_when(
        str_sub(site, 1 , 1) == "H" ~ "high elevation",
        str_sub(site, 1 , 1) == "L" ~ "low elevation",
        str_sub(site, 1 , 1) == "S" ~ "mid elevation"
    )) %>%
    select(population, site, latitude_dec, longitude_dec, elevation_m)

# phila
temp <- sites_phila %>%
    st_as_sf(coords = c("longitude_dec", "latitude_dec"), crs = 4326) %>%
    get_elev_point()

sites_phila <- sites_phila %>%
    mutate(elevation_m = temp$elevation) %>%
    select(population, site, latitude_dec, longitude_dec, elevation_m)

# bind rows
sites <- bind_rows(mutate(sites_mlbs, gradient = "elevation"), mutate(sites_phila, gradient = "urbanization")) %>%
    select(gradient, everything()) %>%
    mutate(across(c(latitude_dec, longitude_dec, elevation_m), function (x) {round(x, 2)}))

# Compute pairwise geo distance
coords_sf <- st_as_sf(sites, coords = c("longitude_dec", "latitude_dec"), crs = 4326)  # WGS84
m <- st_distance(coords_sf)
colnames(m) <- sites$site

sites_dist <- as_tibble(m) %>%
    mutate(site1 = sites$site) %>%
    pivot_longer(-site1, names_to = "site2", values_to = "dist_geo_m") %>%
    left_join(select(sites, site1 = site, gradient1 = gradient)) %>%
    left_join(select(sites, site2 = site, gradient2 = gradient)) %>%
    select(site1, site2, dist_geo_m) %>%
    mutate(dist_geo_km = dist_geo_m / 1000)

write_csv(sites, paste0(folder_phenotypes, "sites/sites.csv"))
write_csv(sites_dist, paste0(folder_phenotypes, "sites/sites_dist.csv"))

# 1. Extract the climatology data from Daymet ----
list_dm <- rep(list(NA), nrow(sites))
for (i in 1:nrow(sites)) {
    dm <- download_daymet(site = sites$site[i],
                          lat = sites$latitude_dec[i],
                          lon = sites$longitude_dec[i],
                          start = 2022,
                          end = 2022,
                          internal = TRUE)
    list_dm[[i]] <- as_tibble(dm$data) %>% mutate(site = sites$site[i], gradient = sites$gradient[i])

}

dml <- bind_rows(list_dm) %>%
    clean_names() %>%
    mutate(ydate = strptime(paste("2022", yday), format="%Y %j")) %>%
    mutate(ymonth = month(ydate)) %>%
    left_join(sites) %>%
    select(gradient, population, site, everything())

write_csv(dml, paste0(folder_phenotypes, "sites/dml.csv"))


# 2. Resample the temperature difference between paired sites ----
tb <- crossing(gradient = c("elevation", "urbanization"), variable = c("tmax_deg_c", "tmin_deg_c"), yday = 1:365)

compute_diff <- function(dml, gra, y, variable) {
    if (gra == "elevation") {
        var1 <- dml %>% filter(gradient == gra, population == "low elevation", yday == y) %>% pull({{variable}})
        var2 <- dml %>% filter(gradient == gra, population == "high elevation", yday == y) %>% pull({{variable}})
    } else if (gra == "urbanization") {
        var1 <- dml %>% filter(gradient == gra, population == "urban", yday == y) %>% pull({{variable}})
        var2 <- dml %>% filter(gradient == gra, population == "suburban", yday == y) %>% pull({{variable}})
    }
    sample_var1 <- sample(var1, n_resample, replace = T)
    sample_var2 <- sample(var2, n_resample, replace = T)
    return(tibble(resample = 1:n_resample,
            sample_var1 = sample_var1,
            sample_var2 = sample_var2,
            diff_var = sample_var1 - sample_var2
        ))
}

set.seed(1)
n_resample <- 100

tbs <- tb %>%
    rowwise() %>%
    mutate(samples = list(compute_diff(dml, gradient, yday, variable)))

diff_vars <- tbs %>%
    ungroup() %>%
    unnest(cols = samples)

write_csv(diff_vars, paste0(folder_phenotypes, "sites/diff_vars.csv"))

# 3. Shade for month ----
tb_month <- dml %>%
    distinct(yday, ydate, ymonth) %>%
    group_by(ymonth) %>%
    filter(ydate == max(ydate) | ydate == min(ydate)) %>%
    mutate(temp = c("start", "end")) %>%
    select(-ydate) %>%
    pivot_wider(names_from = temp, values_from = yday) %>%
    mutate(ymonth = factor(ymonth)) %>%
    ungroup()

write_csv(tb_month, paste0(folder_phenotypes, "sites/tb_month.csv"))
