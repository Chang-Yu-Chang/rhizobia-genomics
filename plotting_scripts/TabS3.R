#' This script compares the traits in pairs of populations inoculated to M. sativa

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

pairs_perm <- read_csv(paste0(folder_data, "phenotypes/plants/sativa/pairs_perm.csv"))
pairs_cont_perm <- read_csv(paste0(folder_data, "phenotypes/plants/sativa_all/pairs_cont_perm.csv"))
pairs_cata_perm <- read_csv(paste0(folder_data, "phenotypes/plants/sativa_all/pairs_cata_perm.csv"))

tb <- bind_rows(pairs_cont_perm, pairs_cata_perm) %>%
    mutate(trait_type = factor(trait_type, c("shoot", "nodule", "leaf", "root"))) %>%
    arrange(gradient, trait_type) %>%
    mutate(ii = rep(1:(n()/2), each = 2))

clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_replace("exp_id", "strain") %>%
        str_replace("exp_labgroup", "labgroup") %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_remove("glmer\\(|lmer\\(")
}

ft <- tb %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab, ii) %>%
    # Clean the table
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        #Trait = factor(Trait, traits$trait_pre),
        P = ifelse(str_detect(Term, "Intercept|sd__"), "", P),
        P = str_replace(P, "p<", "<"),
        P = str_replace(P, "p=", "")
    ) %>%
    select(-ii) %>%
    #arrange(Gradient, Trait) %>%
    flextable() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Type", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Type", "Trait", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Type", "Trait", "Term", "P"), align = "center", part = "all") %>%
    autofit() %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS3.html"), res = 300)
save_as_image(ft, path = here::here("plots/TabS3.png"), res = 300)
