#' This scripts plots the factorial design of the plant experiment, and the plot the manually measured phenotypes

library(tidyverse)
library(cowplot)
library(broom)
library(janitor)
library(ggsci)
library(waffle) #remotes::install_github("hrbrmstr/waffle")
library(lme4) # for linear mixed-effect models
library(car) # Companion to Applied Regression
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)

# 1. Compare H vs. M vs. L vs.  plant fitness using rhizobia strains as environment ----
p <- treatments %>%
    filter(!is.na(dry_weight_mg)) %>%
    filter(rhizobia %in% c("H3M1R1", "L2M2R1")) %>%
    filter(plant_site != "S") %>%
    # mutate(plant_site = case_when(
    #     plant_site == "H" ~ "High",
    #     plant_site == "S" ~ "Mid",
    #     plant_site == "L" ~ "Low",
    # )) %>%
    ggplot(aes(x = plant_site, y = dry_weight_mg, fill = rhizobia_site)) +
    geom_rect(aes(fill = plant_site), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.03) +
    geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, outlier.size = 0, alpha = .7) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), shape = 21, stroke = .5, size = 2) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    facet_grid(~plant_site, scales = "free_x") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "top"
    ) +
    guides(fill = guide_legend(title = "rhizobia site")) +
    labs(x = "plant site", y = "above-ground weight (mg)")
ggsave(paste0(folder_data, "temp/11a-01-biomass_match.png"), p, width = 4, height = 4)

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

# 9. four traits ----
p <- treatments %>%
    drop_na(rhizobia) %>%
    select(id, site = rhizobia_site, dry_weight_mg, nodule_number, root_weight_mg, nodule_weight_mg) %>%
    pivot_longer(cols = c(dry_weight_mg, nodule_number, root_weight_mg, nodule_weight_mg), names_to = "trait") %>%
    ggplot(aes(x = site, y = value, fill = site)) +
    geom_boxplot(alpha = .6, outlier.size = 0, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y", nrow = 1) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        strip.background = element_rect(color = NA, fill = NA)
    ) +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/11a-09-measurements_site.png"), p, width = 8, height = 4)


# 10. for trait in different y scale ----
tt <- treatments %>%
    drop_na(rhizobia) %>%
    select(id, site = rhizobia_site, dry_weight_mg, nodule_number, root_weight_mg, nodule_weight_mg)

p1 <- tt %>%
    drop_na(dry_weight_mg) %>%
    ggplot(aes(x = site, y = dry_weight_mg, fill = site)) +
    geom_boxplot(alpha = .6, outlier.shape = NA, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text.y.left = element_blank(),
        legend.position = "none"
    ) +
    guides() +
    labs(x = "", y = "above-ground dry weight (mg)")

p2 <- tt %>%
    ggplot(aes(x = site, y = root_weight_mg, fill = site)) +
    geom_boxplot(alpha = .6, outlier.shape = NA, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text.x = element_blank(),
        legend.position = "top"
    ) +
    guides() +
    labs(x = "", y = "root dry weight (mg)")

p3 <- tt %>%
    ggplot(aes(x = site, y = nodule_number, fill = site)) +
    geom_boxplot(alpha = .6, outlier.shape = NA, color = "black") +
    geom_jitter(shape = 21, width = 0.2, size = 2, stroke = 1) +
    scale_fill_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text.x = element_blank(),
        legend.position = "none"
    ) +
    guides() +
    labs(x = "", y = "# of nodules")

p <- plot_grid(p1, p2, p3, nrow = 1, axis = "tbrl", align = "hv")

ggsave(paste0(folder_data, "temp/11a-10-measurements_site_trait.png"), p, width = 8, height = 6)







