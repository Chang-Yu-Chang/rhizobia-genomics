#' This script computes the effect sizes of plant experiments
#'
#' 1. Prepare a table for the treatments and continuous response traits
#' 2. Compute effect size for each trait. There are three metrics: Cohen's d. Hedge's g, and partial eta squared

library(tidyverse)
library(janitor)
library(lme4) # for linear mixed-effect models
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
library(vcd) # for computing effect sizes of categorical response
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
    nt = c(rep("without nitrogen", 6+7), rep("with nitrogen", 7), rep("without nitrogen", 7)),
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


if (F) {

# Clean the names
ess <- treatments_eff %>%
    unnest(cohensd) %>%
    clean_names() %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "N-",
        nt == "with nitrogen" ~ "N+"
    )) %>%
    mutate(
        pop = factor(pop, c("elevation", "urbanization")),
        plant = factor(plant, c("lupulina", "sativa"))
    ) %>%
    mutate(response = str_remove(response, "mg") %>% str_replace_all("_", " ")) %>%
    arrange(pop, plant, response)
}
if (F) {

background_df <- tibble(
    pop = factor(c("elevation", "urbanization", "elevation", "urbanization", "elevation"), c("elevation", "urbanization")),
    plant = c("lupulina", "lupulina", "sativa", "sativa", "sativa"),
    host_type = c("source", "source", "alternative", "alternative", "alternative"),
    nt = c("N-", "N-", "N-", "N-", "N+")
)
#
list_treatments <- tibble(
    pop = c(rep("VA", 3), rep("PA", 3), rep("VA", 10), rep("PA", 4)),
    plant = c(rep("lupulina", 6), rep("sativa", 14)),
    nt = c(rep("without nitrogen", 6+5), rep("with nitrogen", 5), rep("without nitrogen", 4)),
    response = c(rep(c("nodule_number", "shoot_biomass_mg", "root_biomass_mg"), 2),
                 rep(c("nodule_number", "shoot_height", "longest_petiole_length", "leaf_number", "leaf_color"), 2),
                 "nodule_number", "shoot_height", "leaf_number", "leaf_color")
)

list_treatments <- list_treatments %>%
    rowwise() %>%
    mutate(tb = list(get(paste0(plant, "s")) %>% filter(exp_nitrogen == nt, population == pop))) %>%
    mutate(
        cohensd = list(compute_cohensd(tb, response)),
        hedgesg = list(compute_hedgesg(tb, response)),
        partialetasquared = list(compute_partialetasquared(tb, response))
    )

}
if (F) {

# Example data
set.seed(123)
study1 <- data.frame(group = rep(c("low", "high"), each = 30), value = c(rnorm(30, mean = 50), rnorm(30, mean = 55)), subject = rep(1:30, 2))
study2 <- data.frame(group = rep(c("low", "high"), each = 30), value = c(rnorm(30, mean = 60), rnorm(30, mean = 65)), subject = rep(1:30, 2))

# Fit linear mixed models
model1 <- lmer(value ~ group + (1|subject), data = study1)
model2 <- lmer(value ~ group + (1|subject), data = study2)

# Cohen's d
d1 <- cohens_d(value ~ group, data = study1)
d2 <- cohens_d(value ~ group, data = study2)

# Eta squared
eta1 <- eta_squared(model1)
eta2 <- eta_squared(model2)

# Print results
print(paste("Cohen's d for study 1:", d1$Cohens_d, "with 95% CI:", d1$CI_low, "-", d1$CI_high))
print(paste("Cohen's d for study 2:", d2$Cohens_d, "with 95% CI:", d2$CI_low, "-", d2$CI_high))

print(paste("Eta squared for study 1:", eta1$Eta2_partial, "with 95% CI:", eta1$CI_low, "-", eta1$CI_high))
print(paste("Eta squared for study 2:", eta2$Eta2_partial, "with 95% CI:", eta2$CI_low, "-", eta2$CI_high))



#
anolelimbs <- read_csv("~/Downloads/anole.csv")

anolelimbs %>%
    ggplot(aes(x = Hurricane, y = Femur)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, shape = 21) +
    theme_classic() +
    theme() +
    guides() +
    labs()

mod1 = lm(Femur ~ Hurricane + SVL, anolelimbs)
Anova(mod1, type=3)
summary(mod1)

# effect sizes
emm.mod1 = emmeans(mod1, specs = "Hurricane")
eff_size(emm.mod1, sigma = sigma(mod1), edf = df.residual(mod1))
cohens_d(Femur ~ Hurricane, data = anolelimbs)

# eta squared
eta_squared(mod1, partial=F)

# partial eta squared
eta_squared(mod1, partial=T)


#
aov(Sepal.Length ~ Species, data = iris) %>%
    parameters::parameters(es_type = c("eta"))

mod <- lm(Sepal.Length ~ Species, data = iris)
Anova(mod)
eta_squared(mod, partial = T)




cohens_d(mpg ~ am, data = mtcars)
cohens_d(mpg ~ am, data = mtcars, pooled_sd = F)
hedges_g(mpg ~ am, data = mtcars)

mean(mtcars$mpg[mtcars$am == 0])
mean(mtcars$mpg[mtcars$am == 1])

var(mtcars$mpg[mtcars$am == 0])
var(mtcars$mpg[mtcars$am == 1])



model <- lm(mpg ~ cyl * am, data = mtcars)
summary(model)
datawizard::standardize(model) %>% summary()
parameters::standardize_parameters(model)


model <- glm(am ~ cyl + hp,
             family = "binomial",
             data = mtcars
)

parameters::standardize_parameters(model, exponentiate = TRUE)

options(contrasts = c("contr.sum", "contr.poly"))

data("ChickWeight")
# keep only complete cases and convert `Time` to a factor
ChickWeight <- subset(ChickWeight, ave(weight, Chick, FUN = length) == 12)
ChickWeight$Time <- factor(ChickWeight$Time)

model <- aov(weight ~ Diet * Time + Error(Chick / Time),
             data = ChickWeight
)
summary(model)
eta_squared(model, partial = TRUE)












}
if (F) {



# Compute effect sizes of categorical responses/traits ----

# Clean the data set
sativas_va_cat <- sativas_va %>%
    select(population, genome_id, site_group, nitrogen_treatment, starts_with("nodule")) %>%
    drop_na(nodule_shape)

# Function to compute Cramér's V
compute_cramersv <- function(tb, response = "nodule_shape", indices) {
    d <- tb[indices, ] # allows boot to select sample
    table_data <- table(d$site_group, d[[response]])
    cramers_v <- assocstats(table_data)$cramer
    return(cramers_v)
}
compute_cramersv_boot <- function (tb, response) {
    boot_results <- boot(tb, function(tb, indices) compute_cramersv(tb, indices, response = response), R = 1000)
    bci <- boot.ci(boot_results, type = "perc")
    cv <- tibble(estimate = bci$t0, lower_cl = bci$percent[4], upper_cl = bci$percent[5])
    return(cv)
}

set.seed(123)
list_responses <- tibble(
        response = rep(c("nodule_shape", "nodule_size", "nodule_color"), 2),
        nt = c(rep("without nitrogen", each = 3), rep("with nitrogen", each = 3))
    ) %>%
    rowwise() %>%
    mutate(tb = list(filter(sativas_va_cat, nitrogen_treatment == nt))) %>%
    mutate(cv = list(compute_cramersv_boot(tb, response))) %>%
    unnest(cv)

# 4. Plot Cramer's V ----
ess <- list_responses %>%
    mutate(id = 1:n()) %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(study_name = paste0(" ", response))

p <- ess %>%
    ggplot() +
    geom_hline(yintercept = c(0,1), linetype = 2) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
    geom_point(aes(x = id, y = estimate, shape = nt), size = 3, stroke = 1) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_x_reverse(breaks = 1:6, labels = ess$study_name, expand = c(0,.8), position = "top") +
    scale_y_continuous(limits = c(0,1), expand = c(0,.03)) +
    coord_flip() +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        panel.background = element_rect(color = NA, fill = "grey95"),
        plot.background = element_rect(color = NA, fill = "white"),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA, fill = "white"),
        legend.box.spacing = unit(0, "mm")
    ) +
    guides() +
    labs(y = "Cramer's V", title = "Elevation strains with sativa")

ggsave(paste0(folder_data, "phenotypes_analysis/effectsize/04-cramersv.png"), p, width = 4, height = 3)



}

