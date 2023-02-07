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

# Exploratory plots ----
experiments %>%
    filter(!is.na(DryWeight)) %>%
    #replace_na(list(DryWeight = 0))
    ggplot(aes(x = RhizobiaSite, y = DryWeight, color = Rhizobia)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1)) +
    theme_classic()

names(experiments)

# Compare H vs. L plant fitness using rhizobia strains as environment
p1 <- experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite %in% c("H", "L")) %>%
    #replace_na(list(DryWeight = 0)) %>%
    ggplot(aes(x = Rhizobia, y = DryWeight, color = PlantSite)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = 1) +
    theme_classic() +
    ggtitle("High vs. Low plants")

ggsave(here::here("plots/02-H_L_plants.png"), p1, width = 4, height = 3)

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
    ggplot(aes(x = Rhizobia, y = DryWeight, color = Rhizobia)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1) +
    theme_classic() +
    ggtitle("Mid plants")
p2
ggsave(here::here("plots/02-M_plants.png"), p2, width = 6, height = 4)

p3 <- experiments %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    ggplot(aes(x = Rhizobia, y = DryWeight, color = Rhizobia)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1) +
    facet_wrap(Plant ~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    ggtitle("Mid plants")

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

# The effect of rhizobia strain and source of origin on plant biomass
mod <- experiments %>%
    filter(PlantSite == "S") %>%
    lmer(DryWeight ~ RhizobiaSite + Rhizobia + (1|Waterblock), data = .)
Anova(mod)
summary(mod)
vif.lme(mod)



# The effect of crossing H/L rhizobia and H/L plants on plant biomass
mod <- experiments %>%
    filter(PlantSite %in% c("L", "H")) %>%
    filter(!is.na(DryWeight)) %>%
    lmer(DryWeight ~ RhizobiaSite + PlantSite + (1 + Rhizobia| Waterblock), data = .)
RhizobiaSite * (1|Rhizobia)
Anova(mod)
summary(mod)
vif.lme(mod)






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

