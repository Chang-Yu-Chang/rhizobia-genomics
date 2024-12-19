#'

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

pairs_rn_perm <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_perm.csv"))

clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_replace("exp_id", "strain") %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_remove("glmer\\(|lmer\\(")
}

ft <- pairs_rn_perm %>%
    select(ii, gradient, trait_pre, st, term, statistic, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab, ii) %>%
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
    merge_v(j = c("Gradient", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Trait", "Model", "Term"), align = "center", part = "all") %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", j = c("Term", "Chisq", "P"), i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS7.html"), res = 300)
save_as_image(ft, path = here::here("plots/TabS7.png"), res = 300)
