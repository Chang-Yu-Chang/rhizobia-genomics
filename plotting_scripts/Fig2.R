#' This script plots the growth traits

library(tidyverse)
library(cowplot)
library(ggh4x)      # for nested facets
library(lme4)       # for lmer
library(car)        # for anova
library(emmeans)    # for emmeans

source(here::here("metadata.R"))

# Prepare the data
iso <- read_csv(paste0(folder_data, "output/iso.csv")) %>%
    mutate(contig_species = factor(contig_species, c("S. meliloti", "S. medicae", "S. canadensis", "S. adhaerens")))
gcl_smooth <- read_csv(paste0(folder_phenotypes, 'growth/gcl_smooth.csv')) # Growth traits per well
gts <- read_csv(paste0(folder_phenotypes, 'growth/gts.csv')) # Growth traits per isolate
gtsl <- gts %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, exp_id, r, lag, maxOD) %>%
    mutate(loglag = -log(lag)) %>%
    pivot_longer(-c(temperature, exp_id), names_to = "trait") %>%
    left_join(select(iso, exp_id, genome_id, contig_species)) %>%
    drop_na(value) %>%
    mutate(symb = case_when(
        contig_species %in% c("S. meliloti", "S. medicae") ~ "symbiotic",
        T ~ "non-symbiotic"
    ))

# Panel A growth curve ----
p1 <- gcl_smooth %>% # Compute mean
    left_join(select(iso, exp_id, genome_id, contig_species)) %>%
    group_by(temperature, exp_id, t) %>%
    mutate(mean_abs = mean(abs_fit)) %>%
    ggplot() +
    geom_line(aes(x = t, y = mean_abs, group = exp_id, color = contig_species), linewidth = .3) +
    geom_text(aes(label = temperature), x = 5, y = .45) +
    scale_color_manual(values = species_colors) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    coord_cartesian(clip = "off") +
    facet_wrap2(~temperature, nrow = 1) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95", linewidth = .5),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 10),
        strip.placement = "outside",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(5, "mm"),
        legend.box.margin = margin(0,0,0,0, "mm"),
        legend.text = element_text(size = 10, face = "italic"),
        plot.background = element_blank()
    ) +
    guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 1))) +
    labs(x = "Time (hour)", y = "O.D.[600nm]")


# Panel B growth traits ----
gtwlm <- gtsl %>%
    #filter(value > 0) %>%
    group_by(temperature, contig_species, trait) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value))) %>%
    group_by(temperature, trait) %>%
    mutate(max_mean_value = max(mean_value, na.rm = T)) %>%
    replace_na(list(ci_value = 0)) %>%
    filter(trait %in% c("r", "lag", "maxOD")) %>%
    mutate(trait = factor(case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "loglag" ~ "lag time (-log(hr))",
        trait == "maxOD" ~ "yield (O.D.[600nm])"
    ), c("growth rate (1/hr)", "lag time (hr)", "yield (O.D.[600nm])")))

p2 <- gtsl %>%
    filter(trait %in% c("r", "lag", "maxOD")) %>%
    mutate(trait = factor(case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "loglag" ~ "lag time (-log(hr))",
        trait == "maxOD" ~ "yield (O.D.[600nm])"
    ), c("growth rate (1/hr)", "lag time (hr)", "yield (O.D.[600nm])"))) %>%
    ggplot() +
    # Each strain
    geom_line(aes(x = temperature, y = value, group = exp_id, color = contig_species), alpha = .1) +
    # Mean value
    geom_ribbon(data = gtwlm, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = contig_species, group = contig_species), inherit.aes = FALSE, alpha = 0.2) +
    geom_line(data = gtwlm, aes(x = temperature, y = mean_value, color = contig_species, group = contig_species)) +
    geom_point(data = gtwlm, aes(x = temperature, y = mean_value, color = contig_species, group = contig_species)) +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    facet_wrap2(~trait, scales = "free_y", nrow = 1, strip.position = "left") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95", linewidth = .5),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 10),
        strip.placement = "outside",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(5, "mm"),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.background = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(color = NA), nrow = 2, direction = "vertical")) +
    labs(x = expression(paste("Temperature (", degree, "C)")))


# ----
p <- plot_grid(
    p1, p2,
    ncol = 1, align = "v", axis = "lr",
    labels = c("A", "B"), scale = .95
) + theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig2.png"), p, width = 8, height = 6)


# Stat ----


## at 25 and 30C, non symbiotic strains have higher growth rate  than symbiotic strains
mod <- gtsl %>%
    filter(temperature == "25c", trait == "r") %>%
    lmer(value ~ symb + (1|contig_species), data = .)
Anova(mod, type = 3)

mod <- gtsl %>%
    filter(temperature == "30c", trait == "r") %>%
    lmer(value ~ symb + (1|contig_species), data = .)
Anova(mod, type = 3)

mod <- gtsl %>%
    filter(temperature == "35c", trait == "r") %>%
    lmer(value ~ symb + (1|contig_species), data = .)
Anova(mod, type = 3)

mod <- gtsl %>%
    filter(temperature == "40c", trait == "r") %>%
    lmer(value ~ symb + (1|contig_species), data = .)
Anova(mod, type = 3)

# for S. meliloti, between 25 and 40C
mod <- gtsl %>%
    filter(trait == "r", contig_species == "S. meliloti") %>%
    lm(value ~ temperature, data = .)
Anova(mod, type = 3)
pairs(emmeans(mod, ~temperature))

mod <- gtsl %>%
    filter(trait == "lag", contig_species == "S. meliloti") %>%
    lm(value ~ temperature, data = .)
Anova(mod, type = 3)
pairs(emmeans(mod, ~temperature))

mod <- gtsl %>%
    filter(trait == "maxOD", contig_species == "S. meliloti") %>%
    lm(value ~ temperature, data = .)
Anova(mod, type = 3)
pairs(emmeans(mod, ~temperature))


## among symbiotic strains
mod <- gtsl %>%
    filter(trait == "r", contig_species %in% c("S. meliloti", "S. medicae")) %>%
    lm(value ~ contig_species + temperature, data = .)
Anova(mod, type = 3)
pairs(emmeans(mod, ~contig_species))

mod <- gtsl %>%
    filter(trait == "lag", contig_species %in% c("S. meliloti", "S. medicae")) %>%
    lm(value ~ temperature + contig_species, data = .)
Anova(mod, type = 3)
pairs(emmeans(mod, ~contig_species))

mod <- gtsl %>%
    filter(trait == "maxOD", contig_species %in% c("S. meliloti", "S. medicae")) %>%
    lm(value ~ temperature + contig_species, data = .)
Anova(mod, type = 3)
pairs(emmeans(mod, ~contig_species))


mod <- gtsl %>%
    filter(temperature == "40c", trait == "r", contig_species %in% c("S. meliloti", "S. medicae")) %>%
    lm(value ~ contig_species, data = .)
summary(mod)

mod <- gtsl %>%
    filter(temperature == "40c", trait == "lag", contig_species %in% c("S. meliloti", "S. medicae")) %>%
    lm(value ~ contig_species, data = .)
summary(mod)

mod <- gtsl %>%
    filter(temperature == "40c", trait == "maxOD", contig_species %in% c("S. meliloti", "S. medicae")) %>%
    lm(value ~ contig_species, data = .)
summary(mod)
