#' This script computes the tradeoff between rhizobia growth vs. symbiosis traits
renv::load()
library(tidyverse)
library(janitor)
source(here::here("metadata.R"))


# Clean data ----
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = sm_species)

gts <- read_csv(paste0(folder_data, "phenotypes_analysis/growth/gts.csv")) %>%
    left_join(isolates) %>%
    clean_names() %>%
    mutate(temperature = factor(temperature, paste0(c(25, 30, 35, 40), "c"))) %>%
    arrange(temperature, exp_id) %>%
    mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    arrange(genome_id) %>%
    select(genome_id, site_group, everything())

lupulinas <- read_csv(paste0(folder_data, "phenotypes_analysis/symbiosis/plants.csv")) %>%
    view
    left_join(isolates) %>%
    clean_names() %>%
    mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    filter(genome_id != "control") %>%
    view
    #filter(population == "VA") %>%
    drop_na(nodule_count) %>%
    arrange(genome_id) %>%
    mutate(nitrogen_treatment = "without nitrogen") %>%
    rename(nodule_number = nodule_count) %>%
    mutate(site_group = factor(site_group, c("low elevation", "high elevation"))) %>%
    select(genome_id, site_group, everything())

sativas <- read_csv(paste0(folder_data, "raw/SymbiosisInSoilData_S24.csv")) %>%
    clean_names() %>%
    rename(genome_id = rhizobia_strain, site_group = elevation) %>%
    mutate(genome_id = tolower(genome_id)) %>%
    left_join(select(isolates, exp_id, genome_id, site)) %>%
    left_join(iso) %>%
    mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    filter(genome_id != "control") %>%
    mutate(site_group = factor(site_group, c("low elevation", "high elevation"))) %>%
    arrange(genome_id) %>%
    select(genome_id, site_group, nitrogen_treatment, everything())

# Compute the strain mean
gts_summ <- gts %>%
    select(genome_id, site_group, temperature, r, lag, max_od) %>%
    pivot_wider(names_from = temperature, values_from = c(r,lag,max_od), names_glue = "gts_{temperature}_{.value}")
t1 <- c("nodule_number", "shoot_biomass_mg", "root_biomass_mg")
t2 <- c("nodule_number", "stem_height", "longest_petiole_length", "leaf_number", "leaf_color")
lupulinas_summ <- lupulinas %>%
    group_by(genome_id, site_group) %>%
    summarize(across(all_of(t1), list(mean = mean, sd = sd), .names = "lup_{.col}_{.fn}")) %>%
    select(ends_with("mean"))
sativas_summ <- sativas %>%
    group_by(genome_id, site_group, nitrogen_treatment) %>%
    filter(nitrogen_treatment == "without nitrogen") %>%
    summarize(across(all_of(t2), list(mean = mean, sd = sd), .names = "sat_{.col}_{.fn}")) %>%
    select(ends_with("mean")) %>%
    #pivot_wider(id_cols = c(genome_id, site_group), names_from = nitrogen_treatment, values_from = ends_with("mean")) %>%
    clean_names()

# Join data
summs <- gts_summ %>%
    left_join(lupulinas_summ) %>%
    left_join(sativas_summ)

names(summs)

summs %>%
    pivot_longer(cols = -c(genome_id, site_group))


# gt vs lupulina
p <- summs %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    drop_na(lup_nodule_number_mean) %>%
    select(genome_id, site_group, starts_with("gts"), starts_with("lup")) %>%
    pivot_longer(cols = starts_with("gts"), names_to = "gts_name", values_to = "gts_value") %>%
    pivot_longer(cols = starts_with("lup"), names_to = "lup_name", values_to = "lup_value") %>%
    #filter(str_detect(gts_name, "_r")) %>%
    ggplot(aes(x = gts_value, y = lup_value)) +
    geom_point(shape = 21, size = 3) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(position = "right") +
    facet_grid(lup_name ~ gts_name, scale = "free", switch = "y") +
    theme_linedraw() +
    theme() +
    guides() +
    labs(x = "rhizobia growth trait", y = "lupulina trait")
ggsave(paste0(folder_data, "phenotypes_analysis/tradeoff/01-gt_vs_lup.png"), p, width = 20, height = 5)

# Sativa vs r
p <- summs %>%
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    drop_na(sat_nodule_number_mean) %>%
    select(genome_id, site_group, starts_with("gts"), starts_with("sat")) %>%
    pivot_longer(cols = starts_with("gts"), names_to = "gts_name", values_to = "gts_value") %>%
    pivot_longer(cols = starts_with("sat"), names_to = "sat_name", values_to = "sat_value") %>%
    #filter(str_detect(gts_name, "r")) %>%
    drop_na() %>%
    ggplot(aes(x = gts_value, y = sat_value)) +
    geom_point(shape = 21, size = 3) +
    geom_smooth(method = "lm") +
    scale_y_continuous(position = "right") +
    facet_grid(sat_name ~ gts_name, scale = "free", switch = "y") +
    theme_linedraw() +
    theme() +
    guides() +
    labs(x = "rhizobia growth trait", y = "sativa trait")
ggsave(paste0(folder_data, "phenotypes_analysis/tradeoff/02-gt_vs_sat.png"), p, width = 20, height = 8)

















