#' PERMANOOVA for plant traits

renv::load()
library(tidyverse)
library(flextable)
library(janitor)
library(broom.mixed)
library(vegan) # for permanova
source(here::here("metadata.R"))

plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
set.seed(1)

tb <- tibble(
    pop = rep(c("VA", "PA"), each = 2),
    host = rep(c("lupulina", "sativa"), 2),
    ff = rep("dm ~ site_group", 4)
) %>%
    mutate(
        dat = list(
            # VA lupulina
            filter(plants, population == "VA", exp_id != "control", exp_plant == "lupulina") %>%
                     select(site_group, shoot_biomass_mg, root_biomass_mg, nodule_number) %>% drop_na,
            # VA sativa
            filter(plants, population == "VA", exp_id != "control",  exp_plant == "sativa", exp_nitrogen == "without nitrogen") %>%
                select(site_group, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>% drop_na,
            # PA lupulina
            filter(plants, population == "PA", exp_id != "control", exp_plant == "lupulina") %>%
                select(site_group, shoot_biomass_mg, root_biomass_mg, nodule_number) %>% drop_na,
            # PA sativa
            filter(plants, population == "PA", exp_id != "control", exp_plant == "sativa") %>%
                select(site_group, shoot_height, nodule_number, leaf_number, leaf_color, lateral_root_number, longest_lateral_root_length) %>% drop_na
        )
    ) %>%
    rowwise() %>%
    mutate(
        pcs = list(as_tibble(prcomp(dat[,-1], scale. = TRUE)$x)),
        dm = list(vegdist(select(pcs, starts_with("PC")), method = "euclidean")),
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
    if(is.na(pv)) {
        return(NA)
    } else if (pv > 0.05) {
        return("-")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv <= 0.001) {
        return("***")
    }
}
tb %>%
    select(pop, p.value) %>%
    rowwise %>%
    mutate(xx = list(detect_sig(p.value))) %>%
    unnest(xx)


ft <- tb %>%
    select(Gradient = pop, Host = host, Predictor = term, df = df, SS = SumOfSqs, R2, `F` = statistic, P = p.value) %>%
    rowwise() %>%
    # Clean the table
    mutate(`Signif.` = list(detect_sig(P))) %>%
    mutate(`F` = round(`F`, 2), SS = round(SS, 2), R2 = round(R2, 2), P = edit_p(P)) %>%
    unnest(`Signif.`) %>%
    mutate(Gradient = ifelse(Gradient == "VA", "elevation", "urbanization")) %>%
    mutate(Host = ifelse(Host == "lupulina", "M. lupulina", "M. sativa")) %>%
    arrange(Gradient) %>%
    #
    flextable() %>%
    valign(valign = "top") %>%
    valign(j = 1:2, valign = "center") %>%
    hline(i = c(3,6,9)) %>%
    autofit() %>%
    compose(part = "header", j = "R2", value = as_paragraph("R", as_sup("2"))) %>%
    bg(bg = "grey90", i = seq(1,11,2), j = 3:9) %>%
    merge_v(j = c("Gradient", "Host")) %>%
    align(j = 1:7, align = "center", part = "all") %>%
    style(part = "body", j = "Host", pr_t = fp_text_default(italic = T)) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS5.html"))

