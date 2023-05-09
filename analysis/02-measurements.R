#' This script analyses the hand measured phenotypes: dry weight, nodule count, root weight, nodule weight

library(tidyverse)
library(broom)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # Companion to Applied Regression
# library(lattice)
# library(lsmeans)
# library(glmmADMB)
source(here::here("analysis/00-metadata.R"))

experiments <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(dry_weight = `dry_weight (mg)`) %>%
    clean_names()
nrow(experiments) # 167 plants
experiments %>% tabyl(rhizobia_site, plant_site, show_missing_levels = T)
experiments %>% tabyl(rhizobia, plant_site, show_missing_levels = T)


#
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), unique(experiments$rhizobia))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")

# 1. Compare H vs. M vs. L vs.  plant fitness using rhizobia strains as environment ----
p <- experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(plant_site = factor(plant_site, c("H", "S", "L"))) %>%
    #filter(plant_site %in% c("H", "L")) %>%
    #replace_na(list(dry_weight = 0)) %>%
    ggplot(aes(x = plant_site, y = dry_weight, fill = rhizobia_site, color = plant_site)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic()
ggsave(here::here("plots/02-01-biomass_match.png"), p, width = 5, height = 4)

##
experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(plant_site %in% c("H", "L")) %>%
    #filter(rhizobia_site == "H") %>%
    glm(dry_weight ~ rhizobia + plant_site, data = .) %>%
    tidy()


# 2. Compare M plant fitness using rhizobia statins as environment ----
p <- experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(plant_site == "S") %>%
    ggplot(aes(x = rhizobia, y = dry_weight, fill = rhizobia_site, alpha = rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = 2) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    guides(alpha = "none")
ggsave(here::here("plots/02-02-biomass_rhizobia.png"), p, width = 6, height = 4)

p <- experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(plant_site == "S") %>%
    ggplot(aes(x = rhizobia, y = dry_weight, fill = rhizobia_site, alpha = rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = 2) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    facet_wrap(plant ~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          panel.grid.major.x = element_line(color = grey(0.5, 0.4), linetype = 2),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    guides(alpha = "none")
ggsave(here::here("plots/02-03-bimass_rhizobia_facet.png"), p, width = 10, height = 6)

## Summary statistics
experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(plant_site == "S") %>%
    group_by(rhizobia) %>%
    summarize(dry_weightMean = mean(dry_weight), dry_weightMedian = median(dry_weight))

##
experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(plant_site == "S") %>%
    aov(dry_weight ~ rhizobia, data = .) %>%
    tidy()

temp <- experiments %>%
    filter(!is.na(dry_weight)) %>%
    filter(plant_site == "S") %>%
    filter(!is.na(rhizobia))

pairwise.t.test(temp$dry_weight, temp$rhizobia, p.adjust.method = "bonferroni")


# 4. nodule number ----
p <- experiments %>%
    #filter(!is.na(nodule_number)) %>%
    #filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(plant_site = factor(plant_site, c("H", "S", "L"))) %>%
    ggplot(aes(x = plant_site, y = nodule_number, fill = rhizobia_site, color = plant_site)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic()
ggsave(here::here("plots/02-04-nodule_interaction.png"), p, width = 5, height = 4)

# 5. nodule number and biomass ----
p <- experiments %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = dry_weight), shape = 21, size = 2, stroke = 1) +
    #scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(here::here("plots/02-05-biomass_nodule.png"), p, width = 4, height = 4)

cor.test(experiments$dry_weight, experiments$nodule_number) %>% tidy()


# 6. nodule number and biomass ----
p <- experiments %>%
    drop_na(rhizobia) %>%
    ggplot() +
    geom_smooth(aes(x = nodule_number, y = dry_weight), method = "lm") +
    geom_point(aes(x = nodule_number, y = dry_weight), shape = 21, size = 1.5, stroke = 1) +
    #scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(rhizobia ~.) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(here::here("plots/02-06-biomass_nodule.png"), p, width = 5, height = 4)

experiments %>%
    drop_na(rhizobia) %>%
    nest(.by = rhizobia) %>%
    mutate(result = map(data, ~cor.test(.x$nodule_number, .x$dry_weight)),
           tided = map(result, tidy)) %>%
    unnest(tided)

# 7. biomass / nodule number ----
p <- experiments %>%
    #filter(!is.na(nodule_number)) %>%
    #filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(plant_site = factor(plant_site, c("H", "S", "L"))) %>%
    ggplot(aes(x = plant_site, y = dry_weight/nodule_number, fill = rhizobia_site, color = plant_site)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic()
ggsave(here::here("plots/02-07-biomass_nodule_interaction.png"), p, width = 5, height = 4)


# 8. biomass / nodule number ----
p <- experiments %>%
    drop_na(dry_weight) %>%
    filter(nodule_number != 0) %>%
    ggplot(aes(x = rhizobia, y = dry_weight/nodule_number, fill = rhizobia_site, alpha = rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = -1) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_y_continuous(expand = c(0.1,0)) +
    theme_classic () +
    guides(alpha = "none")
ggsave(here::here("plots/02-08-biomass_nodule_boxplot.png"), p, width = 6, height = 4)


experiments %>%
    mutate(temp = dry_weight/nodule_number) %>%
    pull(temp) %>% range(na.rm = T)



# Check assumptions ----

# Analysis ----

# The effect of crossing H/L rhizobia and H/L plants on plant biomass
experiments %>%
    filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    lmer(dry_weight ~ rhizobia_site * plant_site + (1|plant_site:Plant) + (1|Waterblock), data = .,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) %>%
    Anova(type = 3)


# The effect of rhizobia strains and origins or thermal environments on plant biomass
experiments %>%
    filter(plant_site == "S") %>%
    lmer(dry_weight ~ rhizobia_site + rhizobia + (1|Waterblock), data = .) %>%
    Anova()

# The effect of rhizobia strains and origins or thermal environments on plant biomass
experiments %>%
    filter(plant_site == "S") %>%
    lmer(dry_weight ~ rhizobia_site * rhizobia + (1|rhizobia_site:rhizobia) + (1|Waterblock), data = .) %>%
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
