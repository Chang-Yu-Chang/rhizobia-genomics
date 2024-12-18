#' This script tidy and combine the data from plant inoculation experiments
#' including traits: shoot biomass, root biomass, and nodule count in the lupulina experiments
#' and traits: leaf, root in the sativa experiment

library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

#' Remove this line after these two strains are sequenced
isolates <- isolates %>%
    bind_rows(tibble(exp_id = c("src-1", "bg-1"), site = c("src", "bg"), population = c("suburban", "urban"), gradient = "urbanization", genome_name = NA, genome_id = c("g_src1", "g_bg1"))) %>%
    filter(!genome_id %in% c("g28", "g20"))


# Lupulina experiments ----
sites <- read_csv(paste0(folder_phenotypes, "sites/sites.csv"))
treatments_ttb <- read_csv(paste0(folder_data, "raw/plants/corrected.csv")) %>% clean_names()
treatments_cyc <- read_csv(paste0(folder_data, "raw/plants/treatments_assigned.csv")) %>% clean_names()

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
    mutate(gradient = "urbanization") %>%
    select(-leaf_number, -height_cm) %>%
    select(id, gradient, exp_id, everything())

treatments_cyc <- treatments_cyc %>%
    rename(exp_id = rhizobia, site = rhizobia_site) %>%
    select(-treatment_id, -label, -nodule_weight_mg) %>%
    replace_na(list(exp_id = "control")) %>%
    mutate(site = str_sub(exp_id, 1, 2)) %>%
    mutate(gradient = "elevation") %>%
    select(id, gradient, exp_id, everything())

# Join the tables
lupulinas <- bind_rows(treatments_cyc, treatments_ttb) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        exp_waterblock = waterblock, exp_plantsite = plant_site, exp_plantmaternal = plant,
        shoot_biomass_mg = dry_weight_mg, root_biomass_mg = root_weight_mg
    ) %>%
    # Join the site and isolate information
    left_join(select(sites, gradient, population, site)) %>%
    left_join(isolates) %>%
    mutate(
        across(c(shoot_biomass_mg, root_biomass_mg), function(x) round(x, 2)),
        exp_id = str_replace(exp_id, "control_\\d", "control"),
        gradient = ifelse(exp_id == "control", "control", gradient),
        population = ifelse(exp_id == "control", "control", population),
        genome_id = ifelse(exp_id == "control", "control", genome_id),
        genome_name = ifelse(exp_id == "control", "control", genome_name),
        population = factor(population, c("high elevation", "low elevation", "urban", "suburban", "control")),
        exp_plant = "lupulina",
        exp_nitrogen = "N-"
    ) %>%
    select(-id) %>%
    select(gradient, population, site, exp_id, genome_id, genome_name, exp_plant, exp_nitrogen, everything())

nrow(lupulinas) # 251 plants
lupulinas %>% filter(!is.na(shoot_biomass_mg), shoot_biomass_mg != 0) %>% nrow # 231 plants with non-zero shoot mass
#write_csv(lupulinas, paste0(folder_phenotypes, "plants/lupulinas.csv"))



# Sativa experiments ----
sativas_va <- readxl::read_xlsx(paste0(folder_data, "raw/plants/SymbiosisInSoilData_S24_groups.xlsx")) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        genome_id = rhizobia_strain,
        population = elevation,
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
        gradient = "elevation",
        exp_plant = "sativa"
    ) %>%
    arrange(genome_id, exp_nitrogen) %>%
    select(gradient, population, site, exp_id, genome_id, genome_name, exp_plant, exp_nitrogen, everything(), -notes)

sativas_pa <- read_csv(paste0(folder_data, "raw/plants/BIOL1102_PooledData.csv")) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        exp_id = rhizobia_strain,
        population = location,
        exp_labgroup = lab_group,
        leaf_color = new_leaf_color
    ) %>%
    mutate(exp_id = tolower(exp_id) %>% str_remove(" "), population = tolower(population),) %>%
    #left_join(iso) %>% # We dont have WGS for bg-1 and src-1
    left_join(isolates) %>%
    mutate(
        site = ifelse(exp_id == "control", "control", site),
        genome_name = ifelse(exp_id == "control", "control", genome_name),
        genome_id = ifelse(exp_id == "control", "control", genome_id),
        genome_id = factor(genome_id, c(isolates$genome_id, "control")),
        gradient = "urbanization",
        exp_plant = "sativa",
        exp_nitrogen = "without nitrogen"
    ) %>%
    arrange(genome_id, exp_nitrogen) %>%
    select(gradient, population, site, exp_id, genome_id, genome_name, exp_nitrogen, exp_labgroup, everything(), -notes, -tube_number)


sativas <- bind_rows(sativas_va, sativas_pa) %>%
    mutate(exp_nitrogen = case_when(
        exp_nitrogen == "without nitrogen" ~ "N-",
        exp_nitrogen == "with nitrogen" ~ "N+"
    )) %>%
    mutate(gradient = ifelse(exp_id == "control", "control", gradient)) %>%
    mutate(population = factor(population, c("low elevation", "high elevation", "urban", "suburban", "control")))


nrow(sativas) # 571 plants
#write_csv(sativas, paste0(folder_phenotypes, "plants/sativas.csv"))

# Bind the two tables ----
plants <- bind_rows(lupulinas, sativas) %>%
    select(gradient, population, site, exp_id, genome_id, genome_name, starts_with("exp_"), everything())
write_csv(plants, paste0(folder_phenotypes, "plants/plants.csv"))

nrow(plants) # 822 plants
plants %>%
    group_by(exp_plant, gradient, population) %>%
    count()

# Groups:   exp_plant, gradient, population [10]
# exp_plant gradient     population         n
# <chr>     <chr>        <fct>          <int>
# 1 lupulina  control      control           24
# 2 lupulina  elevation    high elevation    80
# 3 lupulina  elevation    low elevation     79
# 4 lupulina  urbanization urban             34
# 5 lupulina  urbanization suburban          34
# 6 sativa    control      control          105
# 7 sativa    elevation    high elevation   179
# 8 sativa    elevation    low elevation    182
# 9 sativa    urbanization urban             52
# 10 sativa    urbanization suburban          53
