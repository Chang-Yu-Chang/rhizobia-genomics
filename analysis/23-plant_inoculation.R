#' This script tidy up the data of plant inoculation experiments
#' This script analyses the hand measured phenotypes: dry weight, nodule count, root weight, nodule weight

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 1. raw data wrangling  ----
# 1.1 read data
treatments_ttb <- read_csv(paste0(folder_data, "raw/rhizobia/Terrence/corrected.csv")) %>% clean_names()
treatments_cyc <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>% clean_names()
features_cyc <- read_csv(paste0(folder_data, "raw/rhizobia/05-root_architecture/features.csv"), show_col_types = F)

nrow(treatments_ttb) # 84
nrow(treatments_cyc) # 167

# 1.2 tidy up the variables
treatments_ttb <- treatments_ttb %>%
    rename(id = plant, waterblock = block,
           rhizobia = rhizobia_strain, rhizobia_site = location,
           dry_weight_mg = shoot_weight, nodule_number = nodules, root_weight_mg = root_weight) %>%
    # Correct unit
    mutate(dry_weight_mg = dry_weight_mg*10^3, root_weight_mg = root_weight_mg*10^3) %>%
    # Plant unique id
    mutate(id = id + nrow(treatments_cyc)) %>%
    mutate(rhizobia_site = tolower(rhizobia_site)) %>%
    arrange(id) %>%
    clean_names()

features_cyc <- features_cyc %>%
    clean_names() %>%
    mutate(id = str_replace(file_name, ".png", "") %>% as.numeric())
treatments_cyc <- treatments_cyc %>% left_join(features_cyc, by = join_by(id))



# 1.3 merge the data ----
treatments <- treatments_cyc %>% bind_rows(treatments_ttb) %>%
    replace_na(list(rhizobia = "control", rhizobia_site = "control")) %>%
    mutate(rhizobia = str_replace(rhizobia, "control_\\d", "control")) %>%
    mutate(rhizobia_site = case_when(
        rhizobia_site == "H" ~ "high-elevation",
        rhizobia_site == "L" ~ "low-elevation",
        T ~ rhizobia_site
    )) %>%
    mutate(rhizobia_site = factor(rhizobia_site, c("high-elevation", "low-elevation", "urban", "suburban", "control")))

nrow(treatments) # 251 plants
treatments %>% filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>% nrow # 231 plants with non-zero shoot mass

write_csv(treatments, paste0(folder_data, "temp/23-treatments.csv"))

# 2. reshape data ----

# 2.1 long ----
treatments_long <- treatments %>%
    select(id, rhizobia_site, rhizobia, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0) %>%
    pivot_longer(cols = c(-id, -rhizobia, -rhizobia_site), names_to = "trait")
length(unique(treatments_long$id)) # 231 plants
length(unique(treatments_long$rhizobia)) # 14 rhizobia + 1 control
treatments_long %>%
    arrange(rhizobia_site) %>%
    distinct(id, rhizobia_site, rhizobia) %>%
    group_by(rhizobia_site, rhizobia) %>%
    count()

# 2.2 wide ----
treatments_wide <- treatments %>%
    select(id, rhizobia_site, rhizobia, dry_weight_mg, nodule_number, root_weight_mg) %>%
    filter(!is.na(dry_weight_mg), dry_weight_mg != 0)


write_csv(treatments_long, paste0(folder_data, "temp/23-treatments_long.csv"))
write_csv(treatments_wide, paste0(folder_data, "temp/23-treatments_wide.csv"))















