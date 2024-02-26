#' This script plots the trait data of plant inoculation experiments

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(RColorBrewer)
library(ellipse) # for calculating the ellipse
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(broom.mixed) # for tidy up lme4
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
options(contrasts=c("contr.sum", "contr.poly"))

# Read plant data
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
plants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"))
plants_wide <- read_csv(paste0(folder_data, "temp/23-plants_wide.csv"))

site_group_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", 
                 `suburban` = "#0C6291", `urban` = "#BF4342", control = "grey")

# Read growth rate data
isolates_mapping <- read_csv(paste0(folder_data, "temp/00-isolates_mapping.csv"))
gcs <- read_csv(paste0(folder_data, 'temp/21-gcs.csv'))
gc_summs <- read_csv(paste0(folder_data, 'temp/21-gc_summs.csv'))
gc_prms <- read_csv(paste0(folder_data, 'temp/21-gc_prms.csv'))
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'))

isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    unite(trait, trait, temperature) %>%
    #pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD)) %>%
    left_join(isolates) 

# 1. Plot the 30C 
compute_trait_mean <- function (isolates_gc, tra = "r_30c", pop = "VA") {
igcl <- isolates_gc %>%
    filter(trait == tra) %>%
    filter(population == pop)
igcm <- igcl %>%
    group_by(population, site_group, trait) %>%
    summarize(mean = mean(value, na.rm = T), sem = sd(value, na.rm = T) / sqrt(n()))
    return(list(igcl = igcl, igcm = igcm))
}
plot_dots <- function (igcl, igcm) {
    set.seed(1)
    igcl %>%
        ggplot() +
        geom_rect(data = distinct(igcm, site_group), aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = igcm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = igcm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = unique(igcm$trait))
}

t1 <- compute_trait_mean(isolates_gc)
t2 <- compute_trait_mean(isolates_gc, pop = "PA")
p1 <- plot_dots(t1$igcl, t1$igcm)
p2 <- plot_dots(t2$igcl, t2$igcm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/26a-01-r_30c.png"), p, width = 5, height = 3)

# Does rhizobia sites have effect on shoot biomass?
isolates_test <- filter(isolates_gc, trait == "r_30c", population == "VA") %>% rename(r_30c = value)
mod <- lmer(r_30c ~ site_group + (1|site), data = isolates_test)
Anova(mod, type = 3) # no

# Does rhizobia sites have effect on shoot biomass?
isolates_test <- filter(isolates_gc, trait == "r_30c", population == "PA") %>% rename(r_30c = value)
mod <- lmer(r_30c ~ site_group + (1|site),, data = isolates_test)
Anova(mod, type = 3) # yes


# 2. Plot the thermal response


# Does rhizobia sites have effect on thermal response?
isolates_test <- isolates_gc %>%
    separate(trait, into = c("trait", "temperature")) %>%
    filter(trait == "r") %>% rename(r = value) %>%
    filter(population == "VA")

mod <- lmer(r ~ temperature + site_group + temperature * site_group + (1|site) + (1|exp_id), data = isolates_test)
Anova(mod, type = 3) # no

# Does rhizobia sites have effect on thermal biomass?
isolates_test <- filter(isolates_gc, trait == "r_30c", population == "PA") %>% rename(r_30c = value)
mod <- lmer(r_30c ~ site_group + (1|site),, data = isolates_test)
Anova(mod, type = 3) # yes




# 3. Plot the symbiosis traits comparing the two populations
compute_trait_mean2 <- function (plants_long, tra = "dry_weight_mg", pop = "VA") {
    pl <- plants_long %>%
        filter(trait == tra) %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots2 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(data = distinct(plm, site_group), aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "shoot biomass (mg)")
}

t1 <- compute_trait_mean2(plants_long,tra = "dry_weight_mg", pop = "VA")
t2 <- compute_trait_mean2(plants_long,tra = "dry_weight_mg", pop = "PA")
p1 <- plot_dots2(t1$pl, t1$plm)
p2 <- plot_dots2(t2$pl, t2$plm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_data, "temp/26a-03-shoot_bionass.png"), p, width = 5, height = 3)

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "VA")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # no
# mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # no
# mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # no

# Does rhizobia sites have effect on shoot biomass?
plants_test <- filter(plants, population == "PA")
mod <- lmer(dry_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
Anova(mod, type = 3) # yes
# mod <- lmer(nodule_number ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes
# mod <- lmer(root_weight_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes
