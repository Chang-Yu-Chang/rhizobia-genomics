#' This script tidy and combine the data from plant inoculation experiments
#' including traits: shoot biomass, root biomass, and nodule count in the lupulina experiments
#' and traits: leaf, root in the sativa experiment

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

# Lupulina experiments ----
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv"))
treatments_ttb <- read_csv(paste0(folder_data, "raw/plants/corrected.csv")) %>% clean_names()
treatments_cyc <- read_csv(paste0(folder_data, "raw/plants/treatments_assigned.csv")) %>% clean_names()

nrow(treatments_ttb) # 84 plants
nrow(treatments_cyc) # 167 plants

# Clean up
treatments_ttb <- treatments_ttb %>%
    select(-location) %>%
    rename(
        id = plant, waterblock = block,
        exp_id = rhizobia_strain,
        shoot_biomass_g = shoot_weight,
        root_biomass_g = root_weight
    ) %>%
    # Clean the expid name
    mutate(exp_id = str_replace(exp_id, "b_", "_")) %>%
    mutate(exp_id = str_remove(exp_id, "_c\\d")) %>%
    mutate(exp_id = str_replace(exp_id, "p_", "p")) %>%
    mutate(exp_id = str_replace(exp_id, "_", "-")) %>%
    mutate(exp_id = str_replace(exp_id, "control-\\d", "control")) %>%
    # Add the site name
    mutate(site = str_remove(exp_id, "-\\d")) %>%
    mutate(site = str_remove(site, "\\d$")) %>%
    # Plant unique id
    mutate(id = id + nrow(treatments_cyc)) %>%
    mutate(waterblock = paste0("tbb", waterblock)) %>%
    arrange(id) %>%
    clean_names() %>%
    mutate(population = "PA") %>%
    select(-leaf_number, -height_cm) %>%
    select(id, population, exp_id, everything())

treatments_cyc <- treatments_cyc %>%
    rename(
        shoot_biomass_g = ag_biomass,
        root_biomass_g = bg_biomass,
    ) %>%
    mutate(waterblock = paste0("cyc", block)) %>%
    select(-x8, -x9, -block) %>%
    replace_na(list(exp_id = "control")) %>%
    mutate(site = str_sub(exp_id, 1, 2)) %>%
    mutate(population = "VA") %>%
    select(id, population, site, exp_id, waterblock, everything())

# Join the tables
lupulinas <- bind_rows(treatments_cyc, treatments_ttb) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        exp_waterblock = waterblock,
        exp_plantmaternal = plant
    ) %>%
    # Join the site and isolate information
    left_join(select(sites, population, site)) %>%
    left_join(isolates) %>%
    mutate(
        exp_id = str_replace(exp_id, "control_\\d", "control"),
        #gradient = ifelse(exp_id == "control", "control", gradient),
        population = ifelse(exp_id == "control", "control", population),
        genome_id = ifelse(exp_id == "control", "control", genome_id),
        genome_name = ifelse(exp_id == "control", "control", genome_name),
        #population = factor(population, c("high elevation", "low elevation", "urban", "suburban", "control")),
        exp_plant = "lupulina",
        exp_nitrogen = "N-"
    ) %>%
    select(-id) %>%
    select(population, site, exp_id, genome_id, genome_name, exp_plant, exp_nitrogen, everything())

nrow(lupulinas) # 251 plants
lupulinas %>% filter(!is.na(shoot_biomass_g), shoot_biomass_g != 0) %>% nrow # 231 plants with non-zero shoot mass

# Sativa experiments ----
sativas_va <- readxl::read_xlsx(paste0(folder_data, "raw/plants/SymbiosisInSoilData_S24_groups.xlsx")) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        genome_id = rhizobia_strain,
        exp_nitrogen = nitrogen_treatment,
        shoot_height = stem_height,
        exp_labgroup = group_number,
        exp_labsection = lab_section
    ) %>%
    mutate(genome_id = tolower(genome_id)) %>%
    # Join the site and isolate information
    left_join(isolates) %>%
    mutate(
        site = ifelse(genome_id == "control", "control", site),
        exp_id = ifelse(genome_id == "control", "control", exp_id),
        genome_name = ifelse(genome_id == "control", "control", genome_name),
        genome_id = factor(genome_id, c(isolates$genome_id, "control")),
        shoot_height = shoot_height/10, # mm -> cm
        longest_petiole_length = longest_petiole_length/10, # mm -> cm
        population = "VA",
        exp_plant = "sativa"
    ) %>%
    arrange(genome_id, exp_nitrogen) %>%
    select(population, site, exp_id, genome_id, genome_name, exp_plant, exp_nitrogen, everything(), -notes)

sativas_pa <- read_csv(paste0(folder_data, "raw/plants/BIOL1102_PooledData.csv")) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        exp_id = rhizobia_strain,
        exp_labgroup = lab_group,
        leaf_color = new_leaf_color
    ) %>%
    mutate(exp_id = tolower(exp_id) %>% str_remove(" ")) %>%
    #left_join(iso) %>% # We dont have WGS for bg-1 and src-1
    left_join(isolates) %>%
    mutate(
        site = ifelse(exp_id == "control", "control", site),
        genome_name = ifelse(exp_id == "control", "control", genome_name),
        genome_id = ifelse(exp_id == "control", "control", genome_id),
        genome_id = factor(genome_id, c(isolates$genome_id, "control")),
        population = "PA",
        exp_plant = "sativa",
        exp_nitrogen = "without nitrogen"
    ) %>%
    arrange(genome_id, exp_nitrogen) %>%
    select(population, site, exp_id, genome_id, genome_name, exp_nitrogen, exp_labgroup, everything(), -notes, -tube_number)


sativas <- bind_rows(sativas_va, sativas_pa) %>%
    mutate(exp_nitrogen = case_when(
        exp_nitrogen == "without nitrogen" ~ "N-",
        exp_nitrogen == "with nitrogen" ~ "N+"
    ))


nrow(sativas) # 571 plants

# Bind the two tables ----
plants <- bind_rows(lupulinas, sativas) %>%
    select(population, site, exp_id, genome_id, genome_name, starts_with("exp_"), everything())
write_csv(plants, paste0(folder_phenotypes, "plants/plants.csv"))

nrow(plants) # 822 plants
plants %>%
    group_by(exp_plant, population) %>%
    count()

# exp_plant population     n
# <chr>     <chr>      <int>
# 1 lupulina  PA            68
# 2 lupulina  VA           159
# 3 lupulina  control       24
# 4 sativa    PA           157
# 5 sativa    VA           414
