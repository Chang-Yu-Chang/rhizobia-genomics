#' This script plots the reaction norm of thermal adaptation

library(tidyverse)
library(cowplot)
library(ggh4x)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
iso <- read_csv(paste0(folder_data, "output/iso.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

ztrans <- function (x) (x-mean(x, na.rm = T))/sd(x, na.rm = T)

plants %>%
    left_join(iso) %>%
    filter(exp_plant == "lupulina") %>%
    filter(contig_species == "S. medicae") %>%
    group_by(exp_id, genome_id) %>%
    summarize(shoot_biomass_mg = mean(shoot_biomass_mg, na.rm = T), n())
    #mutate(zt = (shoot_biomass_mg - mean(shoot_biomass_mg, na.rm=T) ) / sd(shoot_biomass_mg, na.rm = T) ) %>%
    #mutate(zt = scale(shoot_biomass_mg)) %>%
    mutate(zt = (shoot_biomass_mg - min(shoot_biomass_mg, na.rm=T) ) / (max(shoot_biomass_mg, na.rm = T) - min(shoot_biomass_mg, na.rm = T))) %>%
    summarize(var_sb = var(zt, na.rm = T))

gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

gtw %>%
    filter(temperature == "30c") %>%
    group_by(exp_id) %>%
    summarize(r = mean(r, na.rm = T)) %>%
    mutate(zt = (r-min(r))/(max(r)-min(r))) %>%
    summarize(var_sb = var(zt, na.rm = T))
