#' This script tidy up the data of plant inoculation experiments
#' This script analyses the hand measured phenotypes: dry weight, nodule count, root weight, nodule weight

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

sites <- read_csv(paste0(folder_data, "temp/22-sites.csv"), show_col_types = F)
treatments_ttb <- read_csv(paste0(folder_data, "raw/rhizobia/Terrence/corrected.csv")) %>% clean_names()
treatments_cyc <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv")) %>% clean_names()
features_cyc <- read_csv(paste0(folder_data, "raw/rhizobia/05-root_architecture/features.csv"))

nrow(treatments_ttb) # 84 plants
nrow(treatments_cyc) # 167 plants

# Clean up
treatments_ttb <- treatments_ttb %>%
    select(-location) %>%
    rename(id = plant, waterblock = block,
           exp_id = rhizobia_strain,
           dry_weight_mg = shoot_weight, nodule_number = nodules, root_weight_mg = root_weight) %>%
    # Clean the expID name
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
    #mutate(site_group = tolower(site_group)) %>%
    arrange(id) %>%
    clean_names() %>%
    mutate(population = "PA")

#
features_cyc <- features_cyc %>%
    clean_names() %>%
    mutate(id = str_replace(file_name, ".png", "") %>% as.numeric())
treatments_cyc <- treatments_cyc %>%
    left_join(features_cyc, by = join_by(id)) %>%
    rename(exp_id = rhizobia, site = rhizobia_site) %>%
    replace_na(list(exp_id = "control")) %>%
    mutate(site = str_sub(exp_id, 1, 2)) %>%
    mutate(population = "VA")


# 1.3 merge the data
plants <- bind_rows(treatments_cyc, treatments_ttb) %>%
    mutate(exp_id = str_replace(exp_id, "control_\\d", "control")) %>%
    left_join(select(sites, population, site, site_group)) %>%
    replace_na(list(site_group = "control")) %>%
    mutate(site_group = factor(site_group, c("high elevation", "low elevation", "urban", "suburban", "control")))

nrow(plants) # 251 plants
plants %>% filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>% nrow # 231 plants with non-zero shoot mass

write_csv(plants, paste0(folder_data, "temp/23-plants.csv"))

# 2. Reshape data
# Long
plants_long <- plants %>%
    select(id, population, site_group, site, exp_id, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>%
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
plants_wide <- plants %>%
    select(id, population, site_group, site, exp_id, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0)

write_csv(plants_long, paste0(folder_data, "temp/23-plants_long.csv"))
write_csv(plants_wide, paste0(folder_data, "temp/23-plants_wide.csv"))

