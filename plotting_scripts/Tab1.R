#' This script computes the statistics for LMM models

renv::load()
library(tidyverse)
library(flextable)
#library(knitr)
#library(KableExtra)
library(janitor)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for permanova
source(here::here("metadata.R"))
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>%
    slice(1:32)

set.seed(1)

gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv")) %>%
    clean_names() %>%
    filter(temperature != "40c") %>%
    #left_join(isolates) %>%
    select(exp_id, r, lag, max_od, temperature) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, max_od)) %>%
    drop_na %>%
    left_join(isolates)

# Master table for growth traits
tb1 <- tibble(
    pop = c("VA", "PA"),
    ff = c("r_30c ~ site_group + (1|site)",
           "r_30c ~ site_group + (1|site)")
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gts, population == pop)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3))[2,])
    ) %>%
    unnest(cols = mod_tided)

# Master table for lupulina traits
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))

tb2 <- tibble(
    pop = c("VA", "PA"),
    ff = c("shoot_biomass_mg ~ site_group + (1|genome_id)",
           "shoot_biomass_mg ~ site_group + (1|genome_id)")
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(plants, population == pop) %>%
                       filter(exp_plant == "lupulina", exp_id != "control") %>%
                       drop_na(shoot_biomass_mg)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3))[2,])
    ) %>%
    unnest(cols = mod_tided)

# Master table for sativa traits
tb3 <- tibble(
    pop = c("VA", "PA"),
    ff = c("shoot_height ~ site_group + (1|genome_id)",
           "shoot_height ~ site_group + (1|genome_id)")
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(plants, population == pop) %>%
                       filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen") %>%
                       drop_na(shoot_height)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3))[2,])
    ) %>%
    unnest(cols = mod_tided)


# Flextable
edit_p <- function (pv) {
    if (pv > 0.05) {
        return(as.character(pv))
    } else if (pv < 0.001) {
        return("<0.001")
    }
}

ft <- bind_rows(tb1, tb2, tb3) %>%
    select(Gradient = pop, Formula = ff, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
    rowwise() %>%
    # Clean the table
    mutate(Chisq = round(Chisq, 2), P = round(P, 3)) %>%
    mutate(P_text = edit_p(P)) %>%
    arrange(Gradient) %>%
    flextable() %>%
    #compose(part = "header", j = "Chisq", value = as_paragraph(as_equation("\\chi^2"))) %>%
    valign(valign = "top") %>%
    style(j = "Formula", pr_t = fp_text_default(italic = T)) %>%
    # P value
    bold(~P<0.05, j = "P_text") %>%
    delete_columns("P") %>%
    compose(part = "header", j = "P_text", value = as_paragraph("P")) %>%
    style(part = "header", j = c("df", "P_text"), pr_t = fp_text_default(italic = T)) %>%
    autofit() %>%
    bg(bg = "grey90", i = c(1,3,5), j = 2:6) %>%
    merge_v(j = "Gradient") %>%
    fix_border_issues()
ft

save_as_html(ft, path = here::here("plots/Tab1.html"))


