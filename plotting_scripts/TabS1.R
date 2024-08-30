#' This script computes the statistics for growth traits ~ site_group for each temp

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

gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv")) %>%
    clean_names() %>%
    replace_na(list(max_od = 0)) %>%
    select(exp_id, r, lag, max_od, temperature, well) %>%
    left_join(isolates)


tb <- tibble(
    pop = rep(c("VA", "PA"), each = 12),
    res = rep(rep(c("r", "lag", "yield"), each = 4), 2),
    temp = rep(c("25c", "30c", "35c", "40c"), 6),
    ff = rep(rep(
        c("r ~ site_group + (1|site)",
          "lag ~ site_group + (1|site)",
          "max_od ~ site_group + (1|site)"
        ), each = 4), 2
    )
) %>%
    filter(!(temp == "40c" & res == "lag")) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gtw, population == pop, temperature == temp)),
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
    select(Gradient = pop, Response = res, Temperature = temp, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
    rowwise() %>%
    # Clean the table
    filter(Predictor != "(Intercept)") %>%
    mutate(Chisq = round(Chisq, 2), P = edit_p(P)) %>%
    mutate(`Signif.` = detect_sig(P)) %>%
    mutate(Gradient = ifelse(Gradient == "VA", "elevation", "urbanization")) %>%
    mutate(Response = factor(Response, c("r", "lag", "yield"))) %>%
    mutate(Temperature = factor(Temperature, c("25c", "30c", "35c", "40c"))) %>%
    arrange(Gradient, Response, Temperature) %>%
    #
    flextable() %>%
    valign(valign = "top") %>%
    valign(j = 1:2, valign = "center") %>%
    hline(i = c(4,7,11,15,18,22)) %>%
    autofit() %>%
    bg(bg = "grey90", i = seq(1,22,2), j = 3:8) %>%
    merge_v(j = c("Gradient", "Response")) %>%
    align(j = 1:8, align = "center", part = "all") %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS1.html"))
