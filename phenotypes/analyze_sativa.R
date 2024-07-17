#' This script analyzes the trait data

renv::load()
library(tidyverse)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(broom.mixed) # for tidy up lme4
source(here::here("metadata.R"))

# Clean data
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = sm_species)

ds <- read_csv(paste0(folder_data, "raw/SymbiosisInSoilData_S24.csv")) %>%
    clean_names() %>%
    rename(genome_id = rhizobia_strain) %>%
    mutate(genome_id = tolower(genome_id)) %>%
    left_join(select(isolates, exp_id, genome_id, site)) %>%
    left_join(iso) %>%
    mutate(genome_id = ifelse(elevation == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id)))


# Nodule number ----
p <- ds %>%
    ggplot(aes(x = genome_id, y = nodule_number)) +
    geom_boxplot(aes(fill = species)) +
    geom_jitter(width = 0.2, shape = 21) +
    facet_grid(nitrogen_treatment~site_group, scales = "free", space = "free_x") +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "phenotypes_analysis/sativa/01-nodule_number.png"), p, width = 10, height = 8)

# nodule ~ elevation w/ N
dsw <- ds %>%
    filter(nitrogen_treatment == "with nitrogen") %>%
    filter(genome_id != "control") %>%
    filter(genome_id != "g2")
mod <- lmer(nodule_number ~ elevation + (1|genome_id), data = dsw)
Anova(mod)
# nodule ~ elevation w/o N
dsw <- ds %>%
    filter(nitrogen_treatment == "without nitrogen") %>%
    filter(genome_id != "control") %>%
    filter(genome_id != "g2")
mod <- lmer(nodule_number ~ elevation + (1|genome_id), data = dsw)
Anova(mod)


p <- ds %>%
    filter(genome_id != "control") %>%
    filter(genome_id != "g2") %>%
    ggplot(aes(x = elevation, y = nodule_number)) +
    geom_boxplot() +
    geom_jitter(, shape = 21, width = 0.2) +
    facet_grid(.~nitrogen_treatment, scales = "free", space = "free_x") +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "phenotypes_analysis/sativa/02-nodule_number.png"), p, width = 8, height = 4)

# Quantitative traits ----
# Pivot longer the quantitative traits
dsl <- ds %>%
    filter(genome_id != "control") %>%
    filter(genome_id != "g2") %>%
    select(genome_id, elevation, nitrogen_treatment, nodule_number, stem_height, longest_petiole_length, leaf_number, leaf_color) %>%
    pivot_longer(cols = -c(genome_id, elevation, nitrogen_treatment), names_to = "trait")

p <- dsl %>%
    ggplot(aes(x = elevation, y = value)) +
    geom_boxplot() +
    geom_jitter(shape = 21, width = 0.2) +
    facet_grid(trait ~ nitrogen_treatment, scale = "free") +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "phenotypes_analysis/sativa/03-traits.png"), p, width = 6, height = 12)

# Quantitative traits ~ elevation
compute_lmm <- function (tb, nit, tr) {
    tbf <- tb %>% filter(nitrogen_treatment == nit, trait == tr)
    mod <- lmer(value ~ elevation + (1|genome_id), data = tbf)
    Anova(mod) %>% return()
}

compute_lmm(dsl, "with nitrogen", "nodule_number")
compute_lmm(dsl, "without nitrogen", "nodule_number")

compute_lmm(dsl, "with nitrogen", "leaf_color")
compute_lmm(dsl, "without nitrogen", "leaf_color")

compute_lmm(dsl, "with nitrogen", "leaf_number")
compute_lmm(dsl, "without nitrogen", "leaf_number")

compute_lmm(dsl, "with nitrogen", "longest_petiole_length")
compute_lmm(dsl, "without nitrogen", "longest_petiole_length")


compute_lmm(dsl, "with nitrogen", "stem_height")
compute_lmm(dsl, "without nitrogen", "stem_height")

# Strain level difference ----
ds %>%
    filter(genome_id %in% c("g4", "g13")) %>%
    select(genome_id, elevation, nitrogen_treatment, nodule_number, stem_height, longest_petiole_length, leaf_number, leaf_color) %>%
    pivot_longer(cols = -c(genome_id, elevation, nitrogen_treatment), names_to = "trait") %>%
    ggplot(aes(x = genome_id, y = value)) +
    geom_boxplot(aes(fill = elevation)) +
    geom_jitter(shape = 21, width = 0.2) +
    facet_grid(trait ~ nitrogen_treatment, scale = "free") +
    theme_bw() +
    theme() +
    guides() +
    labs()















