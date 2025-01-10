#'

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

pairs_perm <- read_csv(paste0(folder_phenotypes, "growth/pairs_perm.csv"))

clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_replace("exp_id", "strain") %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_remove("glmer\\(|lmer\\(")
}

ft <- pairs_perm  %>%
    select(ii, gradient, temperature, trait_pre, st, term, statistic, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Temperature = temperature, Model = st, Term = term, Chisq = statistic, P = siglab, ii) %>%
    # Clean the table
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        P = ifelse(str_detect(Term, "Intercept"), "", P)
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Trait", "Temperature", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Temperature", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Trait", "Temperature", "Model", "Term"), align = "center", part = "all") %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", j = c("Term", "Chisq", "P"), i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()


save_as_html(ft, path = here::here("plots/TabS6.html"), res = 300)
save_as_image(ft, path = here::here("plots/TabS6.png"), res = 300)
