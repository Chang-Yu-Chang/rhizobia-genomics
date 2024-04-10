#' This script tidy and combine the data from plant inoculation experiments
#' including traits: shoot biomass, root biomass, and nodule count

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv"))
treatments_ttb <- read_csv(paste0(folder_data, "raw/rhizobia/Terrence/corrected.csv")) %>% clean_names()
treatments_cyc <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv")) %>% clean_names()

nrow(treatments_ttb) # 84 plants
nrow(treatments_cyc) # 167 plants

# Clean up
treatments_ttb <- treatments_ttb %>%
    select(-location) %>%
    rename(id = plant, waterblock = block,
           exp_id = rhizobia_strain,
           dry_weight_mg = shoot_weight, nodule_number = nodules, root_weight_mg = root_weight) %>%
    # Clean the expid name
    mutate(exp_id = str_replace(exp_id, "b_", "_")) %>%
    mutate(exp_id = str_remove(exp_id, "_c\\d")) %>%
    mutate(exp_id = str_replace(exp_id, "p_", "p")) %>%
    mutate(exp_id = str_replace(exp_id, "_", "-")) %>%
    mutate(exp_id = str_replace(exp_id, "control-\\d", "control")) %>%
    # Add the site name
    mutate(site = str_remove(exp_id, "-\\d")) %>%
    mutate(site = str_remove(site, "\\d$")) %>%
    # Correct unit
    mutate(dry_weight_mg = dry_weight_mg*10^3, root_weight_mg = root_weight_mg*10^3) %>%
    # Plant unique id
    mutate(id = id + nrow(treatments_cyc)) %>%
    arrange(id) %>%
    clean_names() %>%
    mutate(population = "PA") %>%
    select(-leaf_number, -height_cm) %>%
    select(id, population, exp_id, everything())

treatments_cyc <- treatments_cyc %>%
    rename(exp_id = rhizobia, site = rhizobia_site) %>%
    select(-treatment_id, -label, -nodule_weight_mg) %>%
    replace_na(list(exp_id = "control")) %>%
    mutate(site = str_sub(exp_id, 1, 2)) %>%
    mutate(population = "VA") %>%
    select(id, population, exp_id, everything())

# Join the tables
plants <- bind_rows(treatments_cyc, treatments_ttb) %>%
    mutate(exp_id = str_replace(exp_id, "control_\\d", "control")) %>%
    left_join(select(sites, population, site, site_group)) %>%
    replace_na(list(site_group = "control")) %>%
    mutate(site_group = factor(site_group, c("high elevation", "low elevation", "urban", "suburban", "control"))) %>%
    select(id, population, exp_id, site_group, site, everything()) %>%
    rename(shoot_biomass_mg = dry_weight_mg, root_biomass_mg = root_weight_mg, nodule_count = nodule_number) %>%
    mutate(across(c(shoot_biomass_mg, root_biomass_mg), function(x) round(x, 2)))

nrow(plants) # 251 plants
plants %>% filter(!is.na(shoot_biomass_mg), shoot_biomass_mg != 0) %>% nrow # 231 plants with non-zero shoot mass

write_csv(plants, paste0(folder_phenotypes, "symbiosis/plants.csv")) # Symbiosis traits per plant


# 2. Reshape data
# Long
plants_long <- plants %>%
    select(id, population, site_group, site, exp_id, shoot_biomass_mg, root_biomass_mg, nodule_count) %>%
    filter(!is.na(shoot_biomass_mg), shoot_biomass_mg != 0) %>%
    pivot_longer(cols = c(-id, -population, -exp_id, -site_group, -site), names_to = "trait")
length(unique(plants_long$id)) # 231 plants
length(unique(plants_long$exp_id)) # 14 rhizobia + 1 control
plants_long %>%
    arrange(site) %>%
    distinct(id, site_group, exp_id) %>%
    group_by(site_group, exp_id) %>%
    dplyr::count()
# Groups:   site_group, exp_id [15]
#    site_group     exp_id      n
#    <fct>          <chr>   <int>
#  1 high elevation H2M3R1     18
#  2 high elevation H3M1R1     41
#  3 high elevation H4M5R1     17
#  4 low elevation  L2M2R1     37
#  5 low elevation  L3M5R1     18
#  6 low elevation  L4M2R2     18
#  7 urban          40th-2      8
#  8 urban          bg-2        8
#  9 urban          pms-2       7
# 10 urban          ppf-2       9
# 11 suburban       crp1-2      8
# 12 suburban       fp1-2       7
# 13 suburban       gp1-1       8
# 14 suburban       src-2       8
# 15 control        control    19

# Wide
# plants_wide <- plants %>%
#     select(id, population, site_group, site, exp_id, shoot_biomass_mg, root_biomass_mg, nodule_count) %>%
#     filter(!is.na(shoot_biomass_mg), shoot_biomass_mg != 0)

write_csv(plants_long, paste0(folder_phenotypes, "symbiosis/plants_long.csv"))
#write_csv(plants_wide, paste0(folder_phenotypes, "symbiosis/plants_wide.csv"))
