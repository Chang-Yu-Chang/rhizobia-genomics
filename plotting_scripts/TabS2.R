#' PERMANOOVA for growth traits

renv::load()
library(tidyverse)
library(flextable)
library(janitor)
library(broom.mixed)
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
set.seed(1)

gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv")) %>%
    clean_names() %>%
    filter(temperature != "40c") %>%
    select(exp_id, r, lag, max_od, temperature) %>%
    left_join(isolates)

tb <- tibble(
    pop = rep(c("VA", "PA"), each = 3),
    temp = rep(c("25c", "30c", "35c"), 2),
    ff = rep("dm ~ site_group", 6)
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gts, population == pop, temperature == temp) %>% select(site_group, starts_with(c("r", "lag", "max"))) %>% drop_na),
        pcs = list(as_tibble(prcomp(dat[,-1], scale. = TRUE)$x)),
        dm = list(vegdist(select(pcs, starts_with("PC")), method = "euclidean")),
        #test = list(mutate(pcs, site_group = dat$site_group))
        mod = list(adonis2(as.formula(ff), data = mutate(pcs, site_group = dat$site_group), permutations = 10000)),
        mod_tided = list(tidy(mod))
    ) %>%
    unnest(cols = mod_tided)

# Flextable
edit_p <- function (pv) {
    if(is.na(pv)) return(NA)
    if (pv < 0.001) {
        return("<0.001")
    } else {
        return(as.character(round(pv, 3)))
    }
}
detect_sig <- function (pv) {
    if(is.na(pv)) return(NA)
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
    select(Gradient = pop, Temperature = temp, Predictor = term, df = df, SS = SumOfSqs, R2, `F` = statistic, P = p.value) %>%
    rowwise() %>%
    # Clean the table
    mutate(`F` = round(`F`, 2), SS = round(SS, 2), R2 = round(R2, 2), P = edit_p(P)) %>%
    mutate(`Signif.` = detect_sig(P)) %>%
    mutate(Gradient = ifelse(Gradient == "VA", "elevation", "urbanization")) %>%
    arrange(Gradient) %>%
    #
    flextable() %>%
    valign(valign = "top") %>%
    valign(j = 1:2, valign = "center") %>%
    hline(i = c(3,6,9,12,15)) %>%
    autofit() %>%
    compose(part = "header", j = "R2", value = as_paragraph("R", as_sup("2"))) %>%
    bg(bg = "grey90", i = seq(1,18,2), j = 3:9) %>%
    merge_v(j = c("Gradient", "Temperature")) %>%
    align(j = 1:7, align = "center", part = "all") %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS2.html"))

