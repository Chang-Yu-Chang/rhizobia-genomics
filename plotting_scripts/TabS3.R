#' This script computes the statistics for growth traits ~ site_group for each temp

renv::load()
library(tidyverse)
library(flextable)
library(janitor)
library(broom)
library(broom.mixed)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
set.seed(1)

gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv")) %>%
    clean_names() %>%
    #filter(temperature != "40c") %>%
    select(exp_id, r, lag, max_od, temperature, well) %>%
    left_join(isolates)

tb <- tibble(
    pop = rep(c("VA", "PA"), each = 3),
    res = rep(c("r", "lag", "yield"), 2),
    ff = rep(
        c("r ~ site_group*temperature + (1|site) + (1|genome_id)",
          "lag ~ site_group*temperature + (1|site) + (1|genome_id)",
          "max_od ~ site_group*temperature + (1|site) + (1|genome_id)"
        ), 2
    )
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gtw, population == pop)),
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

ft <- tb %>%
    select(Gradient = pop, Response = res, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
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
    hline(i = c(3,6,9,12,15)) %>%
    autofit() %>%
    bg(bg = "grey90", i = seq(1,18,2), j = 3:7) %>%
    merge_v(j = c("Gradient", "Response")) %>%
    align(j = 1:7, align = "center", part = "all") %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS3.html"))

