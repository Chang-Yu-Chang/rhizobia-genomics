#' This script plots the growth data

library(tidyverse)
library(cowplot)
library(ggh4x)
library(flextable)
library(ggh4x) # for nested facets
library(broom.mixed) # for tidying the model outputs
library(lme4) # for lmer
library(car) # for anova
library(boot) # for bootstrapping
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))


gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T))
    # mutate(trait = case_when(
    #     trait == "r" ~ "growth rate (1/hr)",
    #     trait == "lag" ~ "lag time (hr)",
    #     trait == "maxOD" ~ "yield [OD]"
    # ))

# 1. Individual traits ----
strips = strip_nested(
    background_y = elem_list_rect(
        color = NA,
        fill = c("white")
    ),
    text_y = elem_list_text(
        angle = 0
    ),
    bleed = T,
    by_layer_y = F,
    clip = "off", size = "variable"
)

p <- gtwl %>%
    drop_na(value) %>%
    ggplot() +
    geom_boxplot(aes(x = population, y = value, fill = population), alpha = 0.5, width = .3, outlier.size = -1) +
    geom_jitter(aes(x = population, y = value, color = population), width = 0, shape = 16, alpha = 0.3) +
    facet_nested_wrap(gradient+trait~., scale = "free", ncol = 1, strip.position = "left", nest_line = element_line(color = "grey90", linewidth = 1), strip = strips) +
    scale_fill_manual(values = population_colors) +
    scale_color_manual(values = population_colors) +
    scale_x_discrete(position = "top") +
    coord_flip() +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        strip.background = element_rect(color = NA)
    ) +
    guides(fill = "none", color = "none") +
    labs()

ggsave(paste0(folder_data, "phenotypes/growth/01-population_pairs.png"), p, width = 6, height = 6)

# Table
gtwl_n <- gtwl %>% filter(temperature == "30c") %>% drop_na(value)

tb <- tibble(
    gradient = rep(c("elevation", "urbanization"), each = 3),
    trait = rep(c("r", "lag", "maxOD"), 2),
    ff = rep(c(
        "value ~ population + (1|exp_id)"
    ), 6)
) %>%
    #filter(!(temperature == "40c" & trait == "lag")) %>%
    mutate(
        dat = map2(gradient, trait, ~filter(gtwl_n, gradient == .x, trait == .y)),
        mod = map2(dat, ff, ~lmer(as.formula(.y), data = .x)),
        mod_tided = map(mod, ~tidy(Anova(.x, type=3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    select(gradient, trait, ff, term, statistic, df, p.value) %>%
    filter(ff != "value ~ population * temperature + (1|exp_id) + (1|well)")

# gtwl_n %>%
#     filter(gradient == "elevation", temperature == "30c") %>%
#     filter(trait == "yield") %>%
#     lmer(value ~ population + (1|exp_id), data = .)

tb_tidied <- tb %>%
    mutate(across(c("statistic", "df", "p.value"), function (x) round(x, 2))) %>%
    select(Gradient = gradient, Trait = trait, Formula = ff, Term = term, everything())

ft <- tb_tidied %>%
    #select(Trait = trait_pre, Model = st, Effect = effect, Term = term, Estimate = t0, `95% CIs` = cis, ` ` = signlab) %>%
    # Clean the table
    #filter(Predictor != "(Intercept)") %>%
    mutate(
        Formula = str_replace(Formula, "value", Trait),
        Trait = factor(Trait, c("r", "lag", "maxOD"))
    ) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    valign(valign = "top") %>%
    merge_v(j = 1:3) %>%
    valign(j = 1:3, valign = "center") %>%
    #align(j = c("Trait", "Effect"), align = "center", part = "all") %>%
    hline(i = seq(2, nrow_part(.), 2)) %>%
    autofit() %>%
    width(j = 2, 1) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "population")) %>%
    #bg(bg = "grey90", j = "Effect", i = ~Effect == "fixed") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    #line_spacing(j = "Model", space = 1.5) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "growth/01-population_pairs_table.png"), res = 300)
