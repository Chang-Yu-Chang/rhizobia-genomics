#' This script computes the effect sizes of plant experiments
#'
#' 1. Prepare a table for the treatments and continuous response traits
#' 2. Compute effect size for each trait. There are three metrics: Cohen's d. Hedge's g, and partial eta squared

library(tidyverse)
library(janitor)
library(lme4) # for linear mixed-effect models
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
source(here::here("metadata.R"))

# Remove control and nonsymbiotic strains
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
plants <- plants %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
lupulinas <- plants %>% filter(exp_plant == "lupulina")
sativas <- plants %>% filter(exp_plant == "sativa")
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = contig_species, exp_id) %>%
    bind_rows(tibble(genome_id = c("g_src1", "g_bg1"), species = NA, exp_id = c("src-1", "bg-1")))

# Prepare the master table of treatments
treatments <- tibble(
    gra = c(rep("elevation", 3), rep("urbanization", 3), rep("elevation", 14), rep("urbanization", 7)),
    plant = c(rep("lupulina", 6), rep("sativa", 21)),
    nt = c(rep("N-", 6+7), rep("N+", 7), rep("N-", 7)),
    response = c(rep(c("nodule_number", "shoot_biomass_mg", "root_biomass_mg"), 2),
                 rep(c("nodule_number", "primary_root_nodule_number", "lateral_root_nodule_number", "shoot_height", "longest_petiole_length", "leaf_number", "leaf_color"), 2),
                 "nodule_number", "shoot_height", "leaf_number", "leaf_color", "primary_root_length", "lateral_root_number", "longest_lateral_root_length")
)

# Trait abbreviations
clean_trait_names <- function (x) {
    str_split(x, pattern = " ")[[1]] %>% str_sub(1,1) %>% paste(collapse = "") %>% toupper() %>% str_pad(width = 4, side = "right", pad = " ")
}
traits <- tibble(
    trait = c("nodule_number", "root_biomass_mg", "shoot_biomass_mg", "lateral_root_nodule_number", "leaf_color", "leaf_number", "longest_petiole_length", "primary_root_nodule_number", "shoot_height", "lateral_root_number", "longest_lateral_root_length", "primary_root_length"),
    response = c("nodule number", "root biomass ", "shoot biomass", "lateral root nodule number", "leaf color", "leaf number", "longest petiole length", "primary root nodule number", "shoot height", "lateral root number", "longest lateral root length", "primary root length"),
    trait_abbr = map_chr(response, clean_trait_names)
)

# Compute effect sizes of continuous responses/traits
subset_plants <- function (x, y, z) {
    #' Subset plants data for each trait
    get(paste0(x, "s")) %>% filter(exp_nitrogen == y, gradient == z)
}
compute_cohensd <- function (tb, response) {
    # Cohen's d
    formu <- paste0(response, " ~ population + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "population")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    return(as_tibble(es))
}
compute_hedgesg <- function (tb, response) {
    # Hedge's g
    formu <- paste0(response, " ~ population + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "population")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    # Compute Hedge's g which corrects for small sample bias
    ns <- as.vector(table(tb$population))
    J <- 1 - (3 / (4 * (sum(ns)) - 9))
    hg <- es %>%
        as_tibble() %>%
        select(contrast, effect.size, lower.CL, upper.CL) %>%
        mutate(effect.size = effect.size * J,
               lower.CL = lower.CL * J,
               upper.CL = upper.CL * J)
    return(hg)
}
compute_partialetasquared <- function (tb, response) {
    # partial eta squared
    formu <- paste0(response, " ~ population + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    eta_squared(mod, partial = T) %>%
        as_tibble %>%
        return()
}

# Compute effect size for each trait
treatments_eff <- treatments %>%
    mutate(
        plants_data = pmap(list(plant, nt, gra), subset_plants),
        cohensd = map2(plants_data, response, compute_cohensd),
        hedgesg = map2(plants_data, response, compute_hedgesg),
        partialetasquared = map2(plants_data, response, compute_partialetasquared),
    ) %>%
    rename(gradient = gra, exp_plant = plant, exp_nitrogen = nt, trait = response) %>%
    left_join(traits)

#
cohensds <- treatments_eff %>%
    select(gradient, exp_plant, exp_nitrogen, trait, trait_abbr, cohensd) %>%
    unnest(cohensd) %>%
    clean_names()

hedgesgs <- treatments_eff %>%
    select(gradient, exp_plant, exp_nitrogen, trait, trait_abbr, hedgesg) %>%
    unnest(hedgesg) %>%
    clean_names()

partialetasquareds <- treatments_eff %>%
    select(gradient, exp_plant, exp_nitrogen, trait, trait_abbr, partialetasquared) %>%
    unnest(partialetasquared) %>%
    clean_names()

write_csv(cohensds, paste0(folder_data, "phenotypes/effectsize/cohensds.csv"))
write_csv(hedgesgs, paste0(folder_data, "phenotypes/effectsize/hedgesgs.csv"))
write_csv(partialetasquareds, paste0(folder_data, "phenotypes/effectsize/partialetasquareds.csv"))
