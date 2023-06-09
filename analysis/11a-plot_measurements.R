#' This scripts plots the factorial design of the plant experiment, and the plot the manually measured phenotypes

library(tidyverse)
library(broom)
library(janitor)
library(waffle) #remotes::install_github("hrbrmstr/waffle")
library(lme4) # for linear mixed-effect models
library(car) # Companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)

# 0. factorial design ----
treatments %>% tabyl(rhizobia_site, plant_site, show_missing_levels = T)
treatments %>% tabyl(rhizobia, plant_site, show_missing_levels = T)

p <- treatments %>%
    mutate(rhizobia_site = ifelse(is.na(rhizobia_site), "control", rhizobia_site),
           rhizobia = ifelse(is.na(rhizobia), "control", rhizobia)) %>%
    #filter(rhizobia != "control") %>%
    mutate(rhizobia_site = factor(rhizobia_site)) %>%
    group_by(rhizobia_site, rhizobia, plant_site) %>%
    count(.drop = F) %>%
    arrange(rhizobia_site, plant_site) %>%
    ggplot() +
    geom_waffle(aes(fill = plant_site, values = n),
                n_rows = 10, color = "white", radius = unit(2, "pt"), size = 0.33,
                na.rm = TRUE, flip = TRUE) +
    #scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Set1")) +
    #scale_y_continuous(labels = function(x) x * 10, expand = c(0,0)) +
    scale_alpha_manual(values = rhizobia_alphas) +
    facet_grid(~rhizobia_site) +
    coord_equal() +
    theme_void() +
    theme(
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/11a-00-factorial_design.png"), p, width = 5, height = 2)


# 2. Compare H vs. M vs. L vs.  plant fitness using rhizobia strains as environment ----
p <- treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(plant_site = factor(plant_site, c("H", "S", "L"))) %>%
    #filter(plant_site %in% c("H", "L")) %>%
    #replace_na(list(dry_weight_mg = 0)) %>%=8
    ggplot(aes(x = plant_site, y = dry_weight_mg, fill = rhizobia_site, color = plant_site)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic()
ggsave(paste0(folder_data, "temp/11a-01-biomass_match.png"), p, width = 5, height = 4)

##
treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(plant_site %in% c("H", "L")) %>%
    #filter(rhizobia_site == "H") %>%
    glm(dry_weight_mg ~ rhizobia + plant_site, data = .) %>%
    tidy()


# 2. Compare M plant fitness using rhizobia statins as environment ----
p <- treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(plant_site == "S") %>%
    ggplot(aes(x = rhizobia, y = dry_weight_mg, fill = rhizobia_site, alpha = rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = 2) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    guides(alpha = "none")
ggsave(paste0(folder_data, "temp/11a-02-biomass_rhizobia.png"), p, width = 6, height = 4)

p <- treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(plant_site == "S") %>%
    ggplot(aes(x = rhizobia, y = dry_weight_mg, fill = rhizobia_site, alpha = rhizobia)) +
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
ggsave(paste0(folder_data, "temp/11a-03-bimass_rhizobia_facet.png"), p, width = 10, height = 6)

## Summary statistics
treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(plant_site == "S") %>%
    group_by(rhizobia) %>%
    summarize(dry_weightMean = mean(dry_weight_mg), dry_weightMedian = median(dry_weight_mg))

##
treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(plant_site == "S") %>%
    aov(dry_weight_mg ~ rhizobia, data = .) %>%
    tidy()

temp <- treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(plant_site == "S") %>%
    filter(!is.na(rhizobia))

pairwise.t.test(temp$dry_weight_mg, temp$rhizobia, p.adjust.method = "bonferroni")


# 4. nodule number ----
p <- treatments %>%
    #filter(!is.na(nodule_number)) %>%
    #filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(plant_site = factor(plant_site, c("H", "S", "L"))) %>%
    ggplot(aes(x = plant_site, y = nodule_number, fill = rhizobia_site, color = plant_site)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic()
ggsave(paste0(folder_data, "temp/11a-04-nodule_interaction.png"), p, width = 5, height = 4)

# 5. nodule number and biomass ----
p <- treatments %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = dry_weight_mg), shape = 21, size = 2, stroke = 1) +
    #scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/11a-05-biomass_nodule.png"), p, width = 4, height = 4)

cor.test(treatments$dry_weight_mg, treatments$nodule_number) %>% tidy()


# 6. nodule number and biomass ----
p <- treatments %>%
    drop_na(rhizobia) %>%
    ggplot() +
    geom_smooth(aes(x = nodule_number, y = dry_weight_mg), method = "lm") +
    geom_point(aes(x = nodule_number, y = dry_weight_mg), shape = 21, size = 1.5, stroke = 1) +
    #scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(rhizobia ~.) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/11a-06-biomass_nodule.png"), p, width = 5, height = 4)

treatments %>%
    drop_na(rhizobia) %>%
    nest(.by = rhizobia) %>%
    mutate(result = map(data, ~cor.test(.x$nodule_number, .x$dry_weight_mg)),
           tided = map(result, tidy)) %>%
    unnest(tided)

# 7. biomass / nodule number ----
p <- treatments %>%
    #filter(!is.na(nodule_number)) %>%
    #filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    mutate(plant_site = factor(plant_site, c("H", "S", "L"))) %>%
    ggplot(aes(x = plant_site, y = dry_weight_mg/nodule_number, fill = rhizobia_site, color = plant_site)) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, lwd = 1, outlier.size = 2, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_color_manual(values = plant_site_colors) +
    theme_classic()
ggsave(paste0(folder_data, "temp/11a-07-biomass_nodule_interaction.png"), p, width = 5, height = 4)


# 8. biomass / nodule number ----
p <- treatments %>%
    drop_na(dry_weight_mg) %>%
    filter(nodule_number != 0) %>%
    ggplot(aes(x = rhizobia, y = dry_weight_mg/nodule_number, fill = rhizobia_site, alpha = rhizobia)) +
    geom_boxplot(lwd = 1, outlier.size = -1) +
    geom_point(position = position_jitter(width = 0.2), shape = 21, stroke = 1, size = 2) +
    scale_alpha_manual(values = rhizobia_alphas) +
    scale_fill_manual(values = rhizobia_site_colors) +
    scale_y_continuous(expand = c(0.1,0)) +
    theme_classic () +
    guides(alpha = "none")
ggsave(paste0(folder_data, "temp/11a-08-biomass_nodule_boxplot.png"), p, width = 6, height = 4)


treatments %>%
    mutate(temp = dry_weight_mg/nodule_number) %>%
    pull(temp) %>% range(na.rm = T)

# 9. root mass ----

names(treatments)














