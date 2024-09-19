#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(broom.mixed)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

set.seed(1)
edit_p <- function (pv) {
    if (pv < 0.001) {
        return("<0.001")
    } else {
        return(as.character(round(pv, 3)))
    }
}
detect_sig <- function (pv) {
    if (pv > 0.05) {
        return("n.s.")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv < 0.001) {
        return("***")
    }
}

# Correct the naming ----
isolates <- isolates %>%
    mutate(gradient = case_when(
        population == "VA" ~ "elevation",
        population == "PA" ~ "urbanization"
    ), .keep = "unused") %>%
    rename(population = site_group)

gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T))

# Compute the trait ~ population X temperature ----
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv")) %>%
    clean_names() %>%
    replace_na(list(max_od = 0)) %>%
    select(exp_id, r, lag, max_od, temperature, well) %>%
    left_join(isolates)

tb_pertemp <- tibble(
    grad = rep(c("elevation", "urbanization"), each = 12),
    res = rep(rep(c("r", "lag", "yield"), each = 4), 2),
    temp = rep(c("25c", "30c", "35c", "40c"), 6),
    ff = rep(rep(
        c("r ~ population + (1|site)",
          "lag ~ population + (1|site)",
          "max_od ~ population + (1|site)"
        ), each = 4), 2
    )
) %>%
    filter(!(temp == "40c" & res == "lag")) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gtw, gradient == grad, temperature == temp)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    rowwise() %>%
    mutate(signif = detect_sig(p.value))

tb_poptemp <- tibble(
    grad = rep(c("elevation", "urbanization"), each = 3),
    res = rep(c("r", "lag", "yield"), 2),
    ff = rep(
        c("r ~ population*temperature + (1|site) + (1|genome_id)",
          "lag ~ population*temperature + (1|site) + (1|genome_id)",
          "max_od ~ population*temperature + (1|site) + (1|genome_id)"
        ), 2
    )
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gtw, gradient == grad)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    rowwise() %>%
    mutate(signif = detect_sig(p.value))

# Reaction norm ----
gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, well, exp_id, r, lag, max_od) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "max_od" ~ "yield [OD]"
    ))

# Compute the mean
gtwlm <- gtwl %>%
    group_by(gradient, temperature, trait, population) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value))) %>%
    group_by(gradient, temperature, trait) %>%
    mutate(max_mean_value = max(mean_value, na.rm = T))
tb_pertemp <- tb_pertemp %>%
    filter(term == "population") %>%
    mutate(trait = case_when(
        res == "r" ~ "growth rate (1/hr)",
        res == "lag" ~ "lag time (hr)",
        res == "yield" ~ "yield [OD]"
    )) %>%
    left_join(distinct(select(gtwlm, grad = gradient, trait, temp = temperature, max_mean_value))) %>%
    ungroup()

tb_poptemp <- tb_poptemp %>%
    filter(term == "population:temperature") %>%
    mutate(trait = case_when(
        res == "r" ~ "growth rate (1/hr)",
        res == "lag" ~ "lag time (hr)",
        res == "yield" ~ "yield [OD]"
    ), edited_p = edit_p(p.value)
    ) %>%
    ungroup()

#


gg <- gtwl %>%
    filter(gradient == "elevation") %>%
    filter(trait == "growth rate (1/hr)")

tt <- filter(gtwlm, gradient == "elevation") %>%
    filter(trait == "growth rate (1/hr)")
tt_30c <- filter(tt, temperature == "30c")

p1 <- gg %>%
    ggplot() +
    # Individual replicates
    geom_line(aes(x = temperature, y = value, group = well, color = population), alpha = 0) +
    geom_point(data = filter(gg, temperature == "30c"), aes(x = temperature, y = value, color = population), position = position_dodge(width = 0.1), shape = 21, size = 3, stroke = 1, alpha = 0.3) +
    # Mean value of 30c
    geom_point(data = tt_30c, aes(x = temperature, y = mean_value, color = population, fill = population), shape = 21, size = 3, stroke = 1, position = position_dodge(width = 0.1)) +
    geom_errorbar(data = tt_30c, aes(x = temperature, ymin = mean_value-ci_value, ymax = mean_value+ci_value, color = population), inherit.aes = FALSE, position = position_dodge(width = 0.1), width = 0.15, linewidth = 1) +
#    geom_segment(data = tt_30c, aes(x = temperature, xend = temperature, y = mean_value-ci_value, yend = mean_value+ci_value, color = population, group = population), inherit.aes = FALSE, position = position_dodge(width = 0.1)) +
    # Mean value
    #geom_ribbon(data = tt, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = population, group = population), inherit.aes = FALSE, alpha = 0.2) +
    #geom_point(data = tt, aes(x = temperature, y = mean_value, color = population, group = population)) +
    #geom_line(data = tt, aes(x = temperature, y = mean_value, color = population, group = population)) +
    # stats per temperature
    #geom_text(data = filter(tb_pertemp, grad == "elevation"), aes(x = temp, y = max_mean_value, label = signif), vjust = -3, size = 3) +
    # stats temp X population
    #geom_text(data = filter(tb_poptemp, grad == "elevation", trait == "growth rate (1/hr)"), aes(label = paste0("P(pop:temp): ", edited_p)), x = Inf, y = Inf, vjust = 1.1, hjust = 1, size = 3) +
    scale_color_manual(values = site_group_colors, name = "population") +
    scale_fill_manual(values = site_group_colors, name = "population") +
    scale_x_discrete(breaks = c("25c", "30c", "35c", "40c"), labels = c(25, 30, 35, 40), expand = c(.05,0)) +
    facet_wrap(~trait, scales = "free_y", nrow = 1, strip.position = "left") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        strip.placement = "outside",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "right",
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs(x = expression(Temperature*degree*C))

p2 <- gg %>%
    ggplot() +
    # Individual replicates
    geom_line(aes(x = temperature, y = value, group = well, color = population), alpha = 0.1) +
    #geom_point(data = filter(gg, temperature == "30c"), aes(x = temperature, y = value, color = population), position = position_dodge(width = 0.1), shape = 21, size = 3, stroke = 1, alpha = 0.3) +
    # Mean value of 30c
    # geom_point(data = tt_30c, aes(x = temperature, y = mean_value, color = population, fill = population), shape = 21, size = 3, stroke = 1, position = position_dodge(width = 0.1)) +
    # geom_errorbar(data = tt_30c, aes(x = temperature, ymin = mean_value-ci_value, ymax = mean_value+ci_value, color = population), inherit.aes = FALSE, position = position_dodge(width = 0.1), width = 0.15, linewidth = 1) +
    # Mean value
    geom_ribbon(data = tt, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = population, group = population), inherit.aes = FALSE, alpha = 0.2) +
    geom_point(data = tt, aes(x = temperature, y = mean_value, color = population, group = population), size = 3) +
    geom_line(data = tt, aes(x = temperature, y = mean_value, color = population, group = population)) +
    # stats per temperature
    #geom_text(data = filter(tb_pertemp, grad == "elevation"), aes(x = temp, y = max_mean_value, label = signif), vjust = -3, size = 3) +
    # stats temp X population
    #geom_text(data = filter(tb_poptemp, grad == "elevation", trait == "growth rate (1/hr)"), aes(label = paste0("P(pop:temp): ", edited_p)), x = Inf, y = Inf, vjust = 1.1, hjust = 1, size = 3) +
    scale_color_manual(values = site_group_colors, name = "population") +
    scale_fill_manual(values = site_group_colors, name = "population") +
    scale_x_discrete(breaks = c("25c", "30c", "35c", "40c"), labels = c(25, 30, 35, 40), expand = c(.05,0)) +
    facet_wrap(~trait, scales = "free_y", nrow = 1, strip.position = "left") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        strip.placement = "outside",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "right",
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides() +
    labs(x = expression(Temperature*degree*C))



ggsave(here::here("fortalk/growth1.png"), p1, width = 5, height = 4)
ggsave(here::here("fortalk/growth2.png"), p2, width = 5, height = 4)
