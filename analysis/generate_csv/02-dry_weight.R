# This script reads the dry weight data

library(tidyverse)
library(broom)
source(here::here("analysis/00-metadata.R"))

biomass <- read_csv(paste0(folder_data, "raw/rhizobia/04-phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(DryWeight = `DryWeight (mg)`)

#

biomass %>%
    filter(!is.na(DryWeight)) %>%
    #replace_na(list(DryWeight = 0))
    ggplot(aes(x = RhizobiaSite, y = DryWeight, color = Rhizobia)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1)) +
    theme_classic()

names(biomass)

# Compare H vs. L plant fitness using rhizobia strains as environment
p1 <- biomass %>%
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
biomass %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite %in% c("H", "L")) %>%
    #filter(RhizobiaSite == "H") %>%
    glm(DryWeight ~ Rhizobia + PlantSite, data = .) %>%
    tidy()

str(biomass)

# Compare M plant fitness using rhizobia statins as environment
p2 <- biomass %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    ggplot(aes(x = Rhizobia, y = DryWeight, color = Rhizobia)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1) +
    theme_classic() +
    ggtitle("Mid plants")
p2
ggsave(here::here("plots/02-M_plants.png"), p2, width = 6, height = 4)

p3 <- biomass %>%
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
biomass %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    group_by(Rhizobia) %>%
    summarize(DryWeightMean = mean(DryWeight), DryWeightMedian = median(DryWeight))

##
biomass %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    aov(DryWeight ~ Rhizobia, data = .) %>%
    tidy()

temp <- biomass %>%
    filter(!is.na(DryWeight)) %>%
    filter(PlantSite == "S") %>%
    filter(!is.na(Rhizobia))

pairwise.t.test(temp$DryWeight, temp$Rhizobia, p.adjust.method = "bonferroni")




























