#' This script computes the statistics for plant traits

renv::load()
library(tidyverse)
library(flextable)
library(janitor)
library(broom)
library(broom.mixed)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
set.seed(1)


# Master table for lupulina traits
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))

# lup
tb1 <- tibble(
    pop = rep(c("VA", "PA"), each = 3),
    host = "M. lupulina",
    res = rep(c("shoot biomass", "root biomass", "nodule number"), 2),
    ff = rep(
        c("shoot_biomass_mg ~ site_group + (1|genome_id) + (1|exp_waterblock)",
          "root_biomass_mg ~ site_group + (1|genome_id) + (1|exp_waterblock)",
          "nodule_number ~ site_group + (1|genome_id) + (1|exp_waterblock)"
        ), 2
    )
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(plants, population == pop) %>%
                       filter(exp_plant == "lupulina", exp_id != "control") %>%
                       drop_na(shoot_biomass_mg)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3)))
    ) %>%
    unnest(cols = mod_tided)


# Master table for sativa traits
tb2 <- tibble(
    pop = rep(c("VA", "PA"), each = 4),
    host = "M. sativa",
    res = rep(c("shoot height", "nodule number", "leaf number", "leaf color"), 2),
    ff = rep(
        c("shoot_height ~ site_group + (1|genome_id)",
          "nodule_number ~ site_group + (1|genome_id)",
          "leaf_number ~ site_group + (1|genome_id)",
          "leaf_color ~ site_group + (1|genome_id)"
        ), 2
    )
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(plants, population == pop) %>%
                       filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen") %>%
                       drop_na(shoot_height)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3)))
    ) %>%
    unnest(cols = mod_tided)

# Flextable
edit_p <- function (pv) {
    if (pv < 0.001) {
        return("<0.001")
    } else {
        return(as.character(round(pv, 3)))
    }
}
detect_sig <- function (pv) {
    if (pv > 0.05) {
        return("-")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv < 0.001) {
        return("***")
    }
}

ft <- bind_rows(tb1, tb2) %>%
    select(Gradient = pop, Host = host, Response = res, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
    rowwise() %>%
    # Clean the table
    filter(Predictor != "(Intercept)") %>%
    mutate(Chisq = round(Chisq, 2), P = edit_p(P)) %>%
    mutate(`Signif.` = detect_sig(P)) %>%
    mutate(Gradient = ifelse(Gradient == "VA", "elevation", "urbanization")) %>%
    arrange(Gradient) %>%
    #
    flextable() %>%
    valign(valign = "top") %>%
    valign(j = 1:2, valign = "center") %>%
    hline(i = c(3,7,10)) %>%
    autofit() %>%
    bg(bg = "grey90", i = seq(1,13,2), j = 3:8) %>%
    merge_v(j = c("Gradient", "Host", "Response")) %>%
    align(j = 1:8, align = "center", part = "all") %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "body", j = "Host", pr_t = fp_text_default(italic = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS4.html"))


