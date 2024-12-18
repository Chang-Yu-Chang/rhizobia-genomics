#' This script compares the traits in pairs of populations inoculated to M. lupulina

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

pairs_perm <- read_csv(paste0(folder_data, "phenotypes/plants/lupulina/pairs_perm.csv"))

clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_replace("exp_id", "strain") %>%
        str_replace("exp_waterblock", "waterblock") %>%
        str_replace("exp_plantmaternal", "plant_maternal") %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_remove("glmer\\(|lmer\\(")
}

ft2 <- pairs_perm %>%
    # Remove the intercept estimates
    mutate(p_value = ifelse(str_detect(term, "Intercept"), "", p_value)) %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab, ii) %>%
    # Clean the table
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        P = ifelse(str_detect(Term, "Intercept|sd__"), "", P)
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Type", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Type", "Trait", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Type", "Trait", "Term"), align = "center", part = "all") %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft2, path = here::here("plots/TabS4.html"), res = 300)
save_as_image(ft2, path = here::here("plots/TabS4.png"), res = 300)
