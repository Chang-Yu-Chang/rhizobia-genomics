# This script reads the dry weight data

library(tidyverse)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # Companion to Applied Regression
# library(lattice)
# library(lsmeans)
# library(glmmADMB)
source(here::here("analysis/00-metadata.R"))

experiments <- read_csv(paste0(folder_data, "raw/rhizobia/04-phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(DryWeight = `DryWeight (mg)`)

experiments %>%
    filter(PlantSite == "S") %>%
    str()

table(experiments$RhizobiaSite, experiments$Rhizobia, experiments$PlantSite)

xtabs(~ RhizobiaSite + Rhizobia + PlantSite, data = experiments)
xtabs(~ RhizobiaSite + PlantSite, data = experiments %>% filter(PlantSite %in% c("H", "L")))

# Exploratory plots ----
experiments %>%
    filter(!is.na(DryWeight)) %>%
    #replace_na(list(DryWeight = 0))
    ggplot(aes(x = RhizobiaSite, y = DryWeight, color = Rhizobia)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1)) +
    theme_classic()


# Colors
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), unique(experiments$Rhizobia))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")

# Compare H vs. M vs. L vs.  plant fitness using rhizobia strains as environment
p1 <- experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(Rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(PlantSite = factor(PlantSite, c("H", "S", "L"))) %>%
    #filter(PlantSite %in% c("H", "L")) %>%
    #replace_na(list(DryWeight = 0)) %>%
    ggplot(aes(x = PlantSite, y = DryWeight, fill = RhizobiaSite, color = PlantSite)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = 1, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic() +
    ggtitle("High vs Mid vs Low plants")
p1
ggsave(here::here("plots/02-H_L_plants.png"), p1, width = 5, height = 4)

##
experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite %in% c("H", "L")) %>%
    #filter(RhizobiaSite == "H") %>%
    glm(DryWeight ~ Rhizobia + PlantSite, data = .) %>%
    tidy()

str(experiments)

# Compare M plant fitness using rhizobia statins as environment
p2 <- experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    ggplot(aes(x = Rhizobia, y = DryWeight, fill = RhizobiaSite, alpha = Rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = 2) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    guides(alpha = "none") +
    ggtitle("Mid plants")
p2
ggsave(here::here("plots/02-M_plants.png"), p2, width = 6, height = 4)

p3 <- experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    ggplot(aes(x = Rhizobia, y = DryWeight, fill = RhizobiaSite, alpha = Rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = 2) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    facet_wrap(Plant ~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          panel.grid.major.x = element_line(color = grey(0.5, 0.4), linetype = 2),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    guides(alpha = "none") +
    ggtitle("Mid plants")
p3
ggsave(here::here("plots/02-M_plants_facet.png"), p3, width = 10, height = 6)

## Summary statistics
experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    group_by(Rhizobia) %>%
    summarize(DryWeightMean = mean(DryWeight), DryWeightMedian = median(DryWeight))

##
experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    aov(DryWeight ~ Rhizobia, data = .) %>%
    tidy()

temp <- experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    filter(!is.na(Rhizobia))

pairwise.t.test(temp$DryWeight, temp$Rhizobia, p.adjust.method = "bonferroni")


# Check assumptions ----

# Analysis ----

# The effect of crossing H/L rhizobia and H/L plants on plant biomass
experiments %>%
    filter(Rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    lmer(DryWeight ~ RhizobiaSite * PlantSite + (1|PlantSite:Plant) + (1|Waterblock), data = .,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) %>%
    Anova(type = 3)


# The effect of rhizobia strains and origins or thermal environments on plant biomass
experiments %>%
    filter(PlantSite == "S") %>%
    lmer(DryWeight ~ RhizobiaSite + Rhizobia + (1|Waterblock), data = .) %>%
    Anova()

# The effect of rhizobia strains and origins or thermal environments on plant biomass
experiments %>%
    filter(PlantSite == "S") %>%
    lmer(DryWeight ~ RhizobiaSite * Rhizobia + (1|RhizobiaSite:Rhizobia) + (1|Waterblock), data = .) %>%
    Anova(type = 3)





# Test ----
# Load script with custom functions to check assumptions of linear models and calculate VIFs
source("~/Dropbox/lab/local-adaptation/data/test/Wood2018/checkAssumptions.R")
experiment1 <- read_csv("~/Dropbox/lab/local-adaptation/data/test/Wood2018/Experiment1.csv", col_types = cols())
experiment2 <- read_csv("~/Dropbox/lab/local-adaptation/data/test/Wood2018/Experiment2.csv", col_types = cols())

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
options(contrasts=c("contr.sum", "contr.poly"))

# THE EFFECT OF NODULES & GALLS ON FITNESS IN CO-INFECTED PLANTS
# Aboveground biomass
mod = lmer(log(AbovegroundMass.g) ~ Nodules + Galls + RootMass.g + scale(Number.of.Nematode.Eggs)
           + Researcher + (1|Block), experiment1)
Anova(mod, type=3)
summary(mod)
check.assumptions(mod) # Assumptions met
vif.lme(mod)
# Excluding plants that received too many nematode eggs does not qualitatively change the answer
mod = lmer(AbovegroundMass.g ~ Nodules + Galls + RootMass.g + scale(Number.of.Nematode.Eggs)
           + Researcher + (1|Block), subset(experiment1, Number.of.Nematode.Eggs<=400))

# Fruit mass
mod = lmer(FruitMass.g ~ Nodules + Galls + RootMass.g + scale(Number.of.Nematode.Eggs)
           + Researcher + (1|Block), experiment1)
Anova(mod, type=3)
summary(mod)
check.assumptions(mod) # Assumptions met
vif.lme(mod)
# Excluding plants that received too many nematode eggs does not qualitatively change the answer
mod = lmer(FruitMass.g ~ Nodules + Galls + RootMass.g + scale(NumEggs) + Researcher + (1|Block),
           subset(experiment1, Number.of.Nematode.Eggs <=400))





if (FALSE) {
#

PlantGrowth %>%
    as_tibble() %>%
    ggplot(aes(x = group, y = weight)) +
    geom_boxplot() +
    geom_jitter() +
    theme_classic()

fit_plant <- aov(weight ~ group, data = PlantGrowth)
summary(fit_plant)
dummy.coef(fit_plant)
coef(fit_plant)

plot(fit_plant, which = 2)

















}

