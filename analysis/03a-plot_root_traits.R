#' This script plots the root traits

library(tidyverse)
library(cowplot)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(DHARMa) # for checking assumptuons of lm4
library(broom.mixed) # for tidy up lme4
library(ggplotify) # for convert base plot into ggplot
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/03-treatments.csv"), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/03-treatments_long.csv"), show_col_types = F)
treatments_scaled <- read_csv(paste0(folder_data, "temp/03-treatments_scaled.csv"), show_col_types = F)
treatments_scaled_long <- read_csv(paste0(folder_data, "temp/03-treatments_scaled_long.csv"), show_col_types = F)

#
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), unique(treatments$rhizobia))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
traits <- c("dry_weight", "nodule_number", "number_of_root_tips", "number_of_branch_points",
            "total_root_length_px", "branching_frequency_per_px", "network_area_px2",
            "average_diameter_px", "median_diameter_px", "maximum_diameter_px",
            "perimeter_px", "volume_px3", "surface_area_px2")

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
#options(contrasts=c("contr.sum", "contr.poly"))


# 1. Traits by rhizobia location ----
p <- treatments_long %>%
    drop_na(value) %>%
    mutate(trait = factor(trait, traits)) %>%
    ggplot(aes(x = rhizobia_site, y = value, color = rhizobia_site)) +
    geom_boxplot(outlier.size = -1) +
    geom_jitter(width = .1, shape = 21) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y") +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/03a-01-trait_site.png"), p, width = 10, height = 8)

##
treatments_scaled2 <- treatments_scaled %>%
    # Excluding the strains that do not nodulate
    filter(!rhizobia %in% c("H2M3R1", "L4M2R2"))
treatments_scaled_long2 <- treatments_scaled_long %>%
    # Excluding the strains that do not nodulate
    filter(!rhizobia %in% c("H2M3R1", "L4M2R2"))

## Check scaled traits values
p <- treatments_long %>%
    mutate(trait = factor(trait, traits)) %>%
    drop_na(value) %>%
    group_by(trait) %>% mutate(n = n()) %>%
    ggplot() +
    geom_histogram(aes(x = value), color = "black", fill = "white") +
    geom_text(aes(x = Inf, y = Inf, label = paste0("n=", n)), hjust = 2, vjust = 2) +
    facet_wrap(~trait, scales = "free_x") +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/03a-01-trait_value.png"), p, width = 15, height = 10)

p <- treatments_scaled_long2 %>%
    mutate(trait = factor(trait, traits)) %>%
    drop_na(value) %>%
    group_by(trait) %>% mutate(n = n()) %>%
    ggplot() +
    geom_histogram(aes(x = value), color = "black", fill = "white") +
    geom_text(aes(x = Inf, y = Inf, label = paste0("n=", n)), hjust = 2, vjust = 2) +
    facet_wrap(~trait, scales = "free_x") +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/03a-01-trait_value_scaled.png"), p, width = 15, height = 10)


## Sanity check: does the biomass differ between plants with rhizobia from the two sites?
mod <- lmer(dry_weight ~ rhizobia_site+ (1|rhizobia), data = treatments_scaled2)
Anova(mod, type = 3) %>% tidy()
#summary(mod)

## GLMM
traits_mod <- treatments_scaled_long2 %>%
    nest(data = -trait) %>%
    mutate(
        mod = map(data, ~ lmer(value ~ rhizobia_site + (1|rhizobia), data = .x)),
        ano = map(mod, ~ Anova(.x, type = 3)),
        tidied = map(ano, tidy),
        sr = map(mod, simulateResiduals)
        #tidied = map(mod, ~ tidy(.x, conf.int = T))
    ) %>%
    unnest(tidied) %>%
    filter(!str_detect(term, "Intercept")) %>%
    select(-data, -mod, -ano)
traits_mod

p_list <- rep(list(NA), 13)
for (i in 1:13) p_list[[i]] <- as.ggplot(~plot(traits_mod$sr[[i]])) + theme(plot.background = element_rect(color = "black"))
p <- plot_grid(plotlist = p_list, nrow = 4, scale = .9, labels = traits, label_x = 0, hjust = 0) + theme(plot.background = element_rect(fill = "white"))
ggsave(paste0(folder_data, "temp/03a-01-check_assumptions.png"), p, width = 30, height = 20)

# 2. biomass vs. nodule number by sites ----
p <- treatments %>%
    drop_na(rhizobia_site) %>%
    ggplot(aes(x = nodule_number, y = dry_weight, color = rhizobia_site)) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/03a-02-trait_site.png"), p, width = 4, height = 3)

##
lmer(dry_weight ~ rhizobia_site + nodule_number + rhizobia_site:nodule_number + (1|rhizobia), data = treatments) %>%
    tidy(conf.int = T) %>%
    filter(effect == "fixed", !str_detect(term, "Intercept")) %>%
    select(-group)

# 3. pairs ----
library(GGally)
p <- treatments %>%
    drop_na(rhizobia_site) %>%
    #select(all_of(traits)) %>%
    ggpairs(
        columns = which(names(treatments) %in% traits),
        # columns = c(9,10,13),
        aes(color = rhizobia_site),
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", shape = 21)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5))
    )

ggsave(paste0(folder_data, "temp/03a-03-trait_pairs.png"), p, width = 15, height = 15)


p <- treatments_scaled %>%
    drop_na(rhizobia_site) %>%
    #select(all_of(traits)) %>%
    ggpairs(
        columns = which(names(treatments) %in% traits),
        # columns = c(9,10,13),
        aes(color = rhizobia_site),
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", shape = 21)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5))
    )

ggsave(paste0(folder_data, "temp/03a-03-trait_pairs_scaled.png"), p, width = 15, height = 15)

##
lmer(dry_weight ~ rhizobia_site + network_area_px2 + rhizobia_site:network_area_px2 + (1|rhizobia), data = treatments) %>%
    tidy(conf.int = T)


# 4. Traits by strain ----
p <- treatments_scaled_long %>%
    drop_na(value) %>%
    mutate(trait = factor(trait, traits)) %>%
    ggplot(aes(x = rhizobia, y = value, color = rhizobia_site)) +
    geom_boxplot() +
    geom_jitter(width = .1, shape = 21) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/03a-04-trait_strain.png"), p, width = 10, height = 8)


## Sanity check: does the biomass differ between plants with rhizobia from the two sites?
mod <- lmer(dry_weight ~ rhizobia + (1|waterblock), data = treatments_scaled)
Anova(mod, type = 3) %>% tidy()
#summary(mod)

## GLMM
traits_mod <- treatments_scaled_long %>%
    nest(data = -trait) %>%
    mutate(
        mod = map(data, ~ lmer(value ~ rhizobia + (1|waterblock), data = .x)),
        ano = map(mod, ~ Anova(.x, type = 3)),
        tidied = map(ano, tidy),
        sr = map(mod, simulateResiduals)
    ) %>%
    unnest(tidied) %>%
    filter(!str_detect(term, "Intercept")) %>%
    select(-data, -mod, -ano)
traits_mod

p_list <- rep(list(NA), 13)
for (i in 1:13) p_list[[i]] <- as.ggplot(~plot(traits_mod$sr[[i]])) + theme(plot.background = element_rect(color = "black"))
p <- plot_grid(plotlist = p_list, nrow = 4, scale = .9, labels = traits, label_x = 0, hjust = 0) + theme(plot.background = element_rect(fill = "white"))
ggsave(paste0(folder_data, "temp/03a-04-check_assumptions.png"), p, width = 30, height = 20)


#options(contrasts=c("contr.sum", "contr.poly"))
## Sanity check
lm(dry_weight ~ rhizobia, data = treatments_scaled) %>%
    summary()

## OLM
traits_mod <- treatments_scaled_long %>%
    mutate(rhizobia = factor(rhizobia, c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2"))) %>%
    nest(data = -trait) %>%
    mutate(
        mod = map(data, ~ lm(value ~ rhizobia, data = .x)),
        tidied = map(mod, tidy)
    ) %>%
    unnest(tidied) %>%
    #filter(!str_detect(term, "Intercept")) %>%
    #filter(p.value < 0.05)
    select(-data, -mod)

## OLM significance
traits_mod %>%
    mutate(significance = case_when(
        p.value < 0.05 ~ "*",
        p.value < 0.01 ~ "**",
        p.value < 0.001 ~ "***",
        TRUE ~ "n.s."
    )) %>%
    select(trait, term, significance) %>%
    group_by(trait) %>%
    pivot_wider(names_from = term, values_from = significance)

#options(contrasts=c("contr.sum", "contr.poly"))

# treatments_scaled_long %>%
#     drop_na() %>%
#     filter(trait == "nodule_number") %>%
#     group_by(rhizobia) %>%
#     summarize(mean(value, na.rm = T))

#contrasts(Salaries$sex)



























