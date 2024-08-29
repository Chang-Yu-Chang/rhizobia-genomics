#' This script computes the effect sizes of plant experiments

renv::load()
library(tidyverse)
library(janitor)
library(ggsci)
library(ggh4x) # for nested facets
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(emmeans) # estimate marginal means
library(effectsize) # companion to Applied Regression
library(vcd) # for computing effect sizes of categorical response
library(boot)
source(here::here("metadata.R"))

#
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))
lupulinas <- read_csv(paste0(folder_phenotypes, "plants/lupulinas.csv"))
sativas <- read_csv(paste0(folder_phenotypes, "plants/sativas.csv"))

# Remove control and nonsymbiontic strains
plants <- plants %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
lupulinas <- lupulinas %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")
sativas <- sativas %>% filter(!genome_id %in% c("g2", "g3", "g15")) %>% filter(genome_id != "control")

if (F) {

# Clean data ----
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    select(genome_id = genome, species = sm_species, exp_id) %>%
    bind_rows(tibble(genome_id = c("g_src1", "g_bg1"), species = NA, exp_id = c("src-1", "bg-1")))

lupulinas <- read_csv(paste0(folder_data, "phenotypes_analysis/symbiosis/plants.csv")) %>%
    left_join(isolates) %>%
    clean_names() %>%
    rename(nodule_number = nodule_count) %>%
    drop_na(nodule_number) %>%
    mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    filter(genome_id != "control") %>%
    # Remove nonsymbiontic strains
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    arrange(genome_id) %>%
    mutate(nitrogen_treatment = "without nitrogen") %>%
    mutate(site_group = factor(site_group, c("low elevation", "high elevation", "urban", "suburban"))) %>%
    select(population, genome_id, site_group, nitrogen_treatment, everything())

sativas_va <- read_csv(paste0(folder_data, "raw/SymbiosisInSoilData_S24.csv")) %>%
    clean_names() %>%
    rename(genome_id = rhizobia_strain, site_group = elevation) %>%
    mutate(genome_id = tolower(genome_id)) %>%
    left_join(select(isolates, exp_id, genome_id, site)) %>%
    left_join(iso) %>%
    rename(shoot_height = stem_height) %>%
    # Remove nonsymbiotic strains
    filter(!genome_id %in% c("g2", "g3", "g15")) %>%
    mutate(population = "VA") %>%
    mutate(genome_id = ifelse(site_group == "control", "control", genome_id)) %>%
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    filter(genome_id != "control") %>%
    arrange(genome_id) %>%
    select(population, genome_id, site_group, nitrogen_treatment, everything(), -notes)

sativas_pa <- read_csv(paste0(folder_data, "raw/BIOL1102_PooledData.csv")) %>%
    clean_names() %>%
    rename(exp_id = rhizobia_strain, site_group = location, leaf_color = new_leaf_color) %>%
    mutate(exp_id = tolower(exp_id) %>% str_remove(" ")) %>%
    mutate(site_group = tolower(site_group)) %>%
    left_join(iso) %>% # We dont have WGS for bg-1 and src-1
    mutate(genome_id = factor(genome_id, c("control", isolates$genome_id))) %>%
    filter(exp_id != "control") %>%
    mutate(population = "PA") %>%
    arrange(genome_id) %>%
    mutate(nitrogen_treatment = "without nitrogen") %>%
    select(population, genome_id, site_group, nitrogen_treatment, everything(), -notes, -tube_number)

sativas <- bind_rows(sativas_va, sativas_pa) %>%
    mutate(site_group = factor(site_group, c("low elevation", "high elevation", "urban", "suburban")))
}


# Compute effect sizes of continuous responses/traits ----
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

# 1. Plot cohen's d ----
ess <- list_treatments %>%
    unnest(cohensd) %>%
    clean_names() %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(pop = case_when(
        pop == "PA" ~ "city",
        pop == "VA" ~ "mountain"
    )) %>%
    mutate(
        study_name = paste0(" ", response),
        pop = factor(pop, c("mountain", "city")),
        plant = factor(plant, c("lupulina", "sativa"))
        ) %>%
    arrange(pop, plant, response) %>%
    mutate(id = 1:n())

background_df <- tibble(pop = factor(rep(c("mountain", "city"), each = 2), c("mountain", "city")), plant = rep(c("sativa", "lupulina"), 2))
strip <- strip_nested(background_y = list(
    NULL, NULL,
    element_rect(fill = alpha(plant_colors[2], 0.2), color = NA),
    element_rect(fill = alpha(plant_colors[1], 0.2), color = NA),
    element_rect(fill = alpha(plant_colors[2], 0.2), color = NA),
    element_rect(fill = alpha(plant_colors[1], 0.2), color = NA)
))

p <- ess %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
    geom_point(aes(x = id, y = effect_size, shape = nt), size = 3, stroke = 1) +
    scale_x_continuous(breaks = 1:20, labels = ess$study_name, expand = c(0,.8)) +
    #scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    #coord_flip() +
    facet_nested(.~pop + plant, scale = "free_x", space = "free_x", strip = strip) +
    #facet_grid(.~pop, scale = "free_x", space = "free_x") +
    #facet_grid(pop ~., scale = "free_y", space = "free_y", switch = "y") +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(c(0,10,0), "mm"),
        strip.text = element_text(size = 10),
        plot.background = element_rect(color = NA, fill = "white"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA, fill = "white")
    ) +
    guides(color = "none") +
    labs(y = "standardized mean difference (Cohen's d)")

ggsave(paste0(folder_phenotypes, "effectsize/01-cohensd.png"), p, width = 8, height = 5)




if (F) {


 # 2. Plot Hedge's g ----
ess <- list_treatments %>%
    unnest(hedgesg) %>%
    clean_names() %>%
    mutate(id = 1:n()) %>%
    mutate(nt = case_when(
        nt == "without nitrogen" ~ "w/o N",
        nt == "with nitrogen" ~ "w/ N"
    )) %>%
    mutate(pop = case_when(
        pop == "PA" ~ "urbanization",
        pop == "VA" ~ "elevation"
    )) %>%
    mutate(study_name = paste0(" ", response))

p <- ess %>%
    ggplot() +
    geom_rect(data = background_df, aes(fill = plant), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = id, xend = id, y = lower_cl, yend = upper_cl), color = "grey10", linewidth = 1) +
    geom_point(aes(x = id, y = effect_size, shape = nt), size = 3, stroke = 1) +
    scale_x_reverse(breaks = 1:20, labels = ess$study_name, expand = c(0,.8), position = "top") +
    scale_y_continuous(limits = c(-3, 3), expand = c(0,.8), breaks = -3:3) +
    scale_shape_manual(values = c(`w/ N` = 16, `w/o N` = 21), labels = c("with nitrogen", "without nitrogen")) +
    scale_fill_manual(values = plant_colors) +
    coord_flip() +
    facet_nested(pop + plant ~., scale = "free_y", space = "free_y", switch = "y", strip = strip) +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = 2, linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.spacing.y = unit(0, "mm"),
        strip.text = element_text(size = 10),
        plot.background = element_rect(color = NA, fill = "white"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white")
    ) +
    guides(color = "none", fill = "none") +
    labs(y = "standardized mean difference (Hedge's g)")

ggsave(paste0(folder_phenotypes, "effectsize/02-hedgesg.png"), p, width = 8, height = 6)
}

# 3. Plot eta squared partial ----
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
ggsave(paste0(folder_phenotypes, "effectsize/03-partialetasquared.png"), p, width = 8, height = 6)






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



