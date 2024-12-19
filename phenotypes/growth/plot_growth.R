#' This script plots the reaction norm of thermal adaptation

library(tidyverse)
library(cowplot)
library(ggh4x)
library(flextable)
source(here::here("metadata.R"))

pairs_anova <- read_csv(paste0(folder_phenotypes, "growth/pairs_anova.csv"))
pairs_perm <- read_csv(paste0(folder_phenotypes, "growth/pairs_perm.csv"))

# Anova table ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: "))
}
ft <- pairs_anova %>%
    select(ii, Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df, P = siglab) %>%
    mutate(
        Model = map2_chr(st, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
    ) %>%
    arrange(Gradient, Trait) %>%
    select(-ii) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Trait", "Term"), align = "center", part = "all") %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "growth/pairs_anova.png"), res = 300)

# Permutation table ----
ft2 <- pairs_perm %>%
    select(ii, Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab) %>%
    mutate(
        Model = map2_chr(st, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft2, path = paste0(folder_phenotypes, "growth/pairs_perm.png"), res = 300)
