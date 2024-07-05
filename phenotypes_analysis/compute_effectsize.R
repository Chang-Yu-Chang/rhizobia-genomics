#' This script computes the effect sizes of plant experiments

renv::load()
library(tidyverse)
library(janitor)
library(ggsci)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
source(here::here("metadata.R"))


# Clean data
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = sm_species)


# gts <- read_csv(paste0(folder_data, "phenotypes_analysis/growth/gts.csv")) %>%
#     left_join(isolates) %>%
#     clean_names() %>%
#     mutate(temperature = factor(temperature, paste0(c(25, 30, 35, 40), "c"))) %>%
#     arrange(temperature, exp_id) %>%
#     mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
#     mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
#     arrange(genome_id) %>%
#     select(genome_id, site_group, everything())

lupulinas <- read_csv(paste0(folder_data, "phenotypes_analysis/symbiosis/plants.csv")) %>%
    left_join(isolates) %>%
    clean_names() %>%
    mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    filter(genome_id != "control") %>%
    filter(population == "VA") %>%
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


# Compute effect sizes
# Cohen's d
compute_cohensd <- function (tb, response) {
    formu <- paste0(response, " ~ site_group + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "site_group")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    return(as_tibble(es))
}
# Hedge's g
compute_hedgesg <- function (tb, response) {
    # Compute Cohen's d
    formu <- paste0(response, " ~ site_group + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    emm.mod <- emmeans(mod, specs = "site_group")
    es <- eff_size(emm.mod, sigma = sigma(mod), edf = df.residual(mod))
    # Compute Hedge's g which corrects for small sample bias
    ns <- as.vector(table(tb$site_group))
    J <- 1 - (3 / (4 * (sum(ns)) - 9))
    hg <- es %>%
        as_tibble() %>%
        select(contrast, effect.size, lower.CL, upper.CL) %>%
        mutate(effect.size = effect.size * J,
               lower.CL = lower.CL * J,
               upper.CL = upper.CL * J)
    return(hg)
}

# partial eta squared
compute_partialetasquared <- function (tb, response) {
    formu <- paste0(response, " ~ site_group + (1|genome_id)")
    mod <- lmer(as.formula(formu), data = tb)
    eta_squared(mod, partial = T) %>%
        as_tibble %>%
        return()
}


#
list_treatments <- tibble(
    plant = c(rep("lupulina", 3), rep("sativa", 10)),
    nt = c(rep("without nitrogen", 3+5), rep("with nitrogen", 5)),
    response = c("nodule_number", "shoot_biomass_mg", "root_biomass_mg",
                 rep(c("nodule_number", "stem_height", "longest_petiole_length", "leaf_number", "leaf_color"), 2))
)
list_treatments <- list_treatments %>%
    rowwise() %>%
    mutate(
        cohensd = list(compute_cohensd(get(paste0(plant, "s")) %>% filter(nitrogen_treatment == nt), response)),
        hedgesg = list(compute_hedgesg(get(paste0(plant, "s")) %>% filter(nitrogen_treatment == nt), response)),
        partialetasquared = list(compute_partialetasquared(get(paste0(plant, "s")) %>% filter(nitrogen_treatment == nt), response))
    )

# 1. Plot cohen's d ----
ess <- list_treatments %>%
    unnest(cohensd) %>%
    clean_names() %>%
    mutate(id = 1:n()) %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(study_name = paste(nt, response, sep = ", "))

p <- ess %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = id, y = effect_size, color = plant), size = 3) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl, color = plant), linewidth = 1) +
    scale_x_reverse(breaks = 1:13, labels = ess$study_name, expand = c(0,.7)) +
    scale_color_aaas() +
    coord_flip() +
    facet_grid(plant ~., scale = "free_y", space = "free_y") +
    theme_light() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 15)
    ) +
    guides(color = "none") +
    labs(y = "standardized mean difference (Cohen's d)")

ggsave(paste0(folder_data, "phenotypes_analysis/effectsize/01-cohensd.png"), p, width = 8, height = 6)

# 2. Plot Hedge's g ----
ess <- list_treatments %>%
    unnest(hedgesg) %>%
    clean_names() %>%
    mutate(id = 1:n()) %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(study_name = paste(nt, response, sep = ", "))

p <- ess %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = id, y = effect_size, color = plant), size = 3) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl, color = plant), linewidth = 1) +
    scale_x_reverse(breaks = 1:13, labels = ess$study_name, expand = c(0,.7)) +
    scale_color_aaas() +
    coord_flip() +
    facet_grid(plant ~., scale = "free_y", space = "free_y") +
    theme_light() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 15)
    ) +
    guides(color = "none") +
    labs(y = "standardized mean difference (Hedge's g)")

ggsave(paste0(folder_data, "phenotypes_analysis/effectsize/02-hedgesg.png"), p, width = 8, height = 6)

# 3. Plot eta squared partial
ess <- list_treatments %>%
    unnest(partialetasquared) %>%
    clean_names() %>%
    mutate(id = 1:n()) %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(study_name = paste(nt, response, sep = ", "))

p <- ess %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = id, y = eta2_partial, color = plant), size = 3) +
    geom_segment(aes(x = id, xend = id, y = ci_low, yend = ci_high, color = plant), linewidth = 1) +
    scale_x_reverse(breaks = 1:13, labels = ess$study_name, expand = c(0,.7)) +
    scale_color_aaas() +
    coord_flip() +
    facet_grid(plant ~., scale = "free_y", space = "free_y") +
    theme_light() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 15)
    ) +
    guides(color = "none") +
    labs(y = "partial eta squared")
ggsave(paste0(folder_data, "phenotypes_analysis/effectsize/03-partialetasquared.png"), p, width = 8, height = 6)






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



