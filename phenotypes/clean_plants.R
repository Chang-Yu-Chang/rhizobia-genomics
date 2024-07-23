#' This script tidy and combine the data from plant inoculation experiments
#' including traits: shoot biomass, root biomass, and nodule count

renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))

# iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
#     select(genome_id = genome, species = sm_species, exp_id) %>%
#     bind_rows(tibble(genome_id = c("g_src1", "g_bg1"), species = NA, exp_id = c("src-1", "bg-1")))

#' Remove this line after these two strains are sequenced
isolates <- isolates %>%
    bind_rows(tibble(exp_id = c("src-1", "bg-1"), site = c("src", "bg"), site_group = c("suburban", "urban"), population = "PA", genome_name = NA, genome_id = c("g_src1", "g_bg1")))


# Lupulina experiments ----
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
lupulinas <- bind_rows(treatments_cyc, treatments_ttb) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        exp_waterblock = waterblock, exp_plantsite = plant_site, exp_plantmaternal = plant,
        shoot_biomass_mg = dry_weight_mg, root_biomass_mg = root_weight_mg
    ) %>%
    # Join the site and isolate information
    left_join(select(sites, population, site, site_group)) %>%
    left_join(isolates) %>%
    mutate(
        across(c(shoot_biomass_mg, root_biomass_mg), function(x) round(x, 2)),
        exp_id = str_replace(exp_id, "control_\\d", "control"),
        site_group = ifelse(exp_id == "control", "control", site_group),
        genome_id = ifelse(exp_id == "control", "control", genome_id),
        genome_name = ifelse(exp_id == "control", "control", genome_name),
        site_group = factor(site_group, c("high elevation", "low elevation", "urban", "suburban", "control")),
        exp_plant = "lupulina",
        exp_nitrogen = "without nitrogen"
    ) %>%
    select(-id) %>%
    select(population, site_group, site, exp_id, genome_id, genome_name, exp_plant, exp_nitrogen, everything())

nrow(lupulinas) # 251 plants
lupulinas %>% filter(!is.na(shoot_biomass_mg), shoot_biomass_mg != 0) %>% nrow # 231 plants with non-zero shoot mass
write_csv(lupulinas, paste0(folder_phenotypes, "plants/lupulinas.csv")) # Symbiosis traits per plant



# Sativa experiments ----
sativas_va <- read_csv(paste0(folder_data, "raw/SymbiosisInSoilData_S24.csv")) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        genome_id = rhizobia_strain, site_group = elevation,
        exp_nitrogen = nitrogen_treatment,
        shoot_height = stem_height
    ) %>%
    mutate(genome_id = tolower(genome_id)) %>%
    # Join the site and isolate information
    left_join(isolates) %>%
    mutate(
        site = ifelse(genome_id == "control", "control", site),
        exp_id = ifelse(genome_id == "control", "control", exp_id),
        genome_name = ifelse(genome_id == "control", "control", genome_name),
        genome_id = factor(genome_id, c(isolates$genome_id, "control")),
        #site_group = factor(site_group, c("high elevation", "low elevation", "urban", "suburban", "control")),
        shoot_height = shoot_height/10,
        population = "VA",
        exp_plant = "sativa"
    ) %>%
    arrange(genome_id, exp_nitrogen) %>%
    select(population, site_group, site, exp_id, genome_id, genome_name, exp_plant, exp_nitrogen, everything(), -notes)

sativas_pa <- read_csv(paste0(folder_data, "raw/BIOL1102_PooledData.csv")) %>%
    # Make sure the variable names are correct
    clean_names() %>%
    rename(
        exp_id = rhizobia_strain, site_group = location,
        exp_labgroup = lab_group,
        leaf_color = new_leaf_color
    ) %>%
    mutate(exp_id = tolower(exp_id) %>% str_remove(" "), site_group = tolower(site_group),) %>%
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
    select(population, site_group, site, exp_id, genome_id, genome_name, exp_nitrogen, exp_labgroup, everything(), -notes, -tube_number)


sativas <- bind_rows(sativas_va, sativas_pa) %>%
    mutate(site_group = factor(site_group, c("low elevation", "high elevation", "urban", "suburban", "control")))


nrow(sativas) # 571 plants
write_csv(sativas, paste0(folder_phenotypes, "plants/sativas.csv"))

# Bind the two tables ----
plants <- bind_rows(lupulinas, sativas)
write_csv(plants, paste0(folder_phenotypes, "plants/plants.csv"))

nrow(plants) # 822 plants
plants %>%
    group_by(exp_plant, population, site_group) %>%
    count()

# Groups:   exp_plant, population, site_group [12]
# exp_plant population site_group         n
# <chr>     <chr>      <fct>          <int>
#     1 lupulina  PA         urban             34
# 2 lupulina  PA         suburban          34
# 3 lupulina  PA         control           16
# 4 lupulina  VA         high elevation    80
# 5 lupulina  VA         low elevation     79
# 6 lupulina  VA         control            8
# 7 sativa    PA         urban             52
# 8 sativa    PA         suburban          53
# 9 sativa    PA         control           52
# 10 sativa    VA         high elevation   179
# 11 sativa    VA         low elevation    182
# 12 sativa    VA         control           53


if (F) {

    #mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    #mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%

    # Remove nonsymbiontic strains
    #filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    #filter(genome_id != "control") %>%
    #drop_na(nodule_number) %>%

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

    write_csv(plants_long, paste0(folder_phenotypes, "plants/plants_long.csv"))
    #write_csv(plants_wide, paste0(folder_phenotypes, "plants/plants_wide.csv"))



}
