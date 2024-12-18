#' This script compares the traits in pairs of populations inoculated to M. sativa

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

# Stat tables
pairs_anova <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_anova.csv"))
pairs_perm <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_perm.csv"))
# For histogram
pairs_perm_obv <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_perm_obv.csv"))
pairs_perm_raw <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_perm_raw.csv"))

clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_replace("exp_id", "strain") %>%
        str_replace("exp_labgroup", "labgroup") %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_remove("glmer\\(|lmer\\(")
}


# Anova table ----
ft2 <- pairs_perm %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab, ii) %>%
    # Clean the table
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        P = ifelse(str_detect(Term, "Intercept|sd__"), "", P),
        P = str_replace(P, "p<", "<"),
        P = str_replace(P, "p=", "")
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
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

save_as_html(ft2, path = here::here("plots/TabS3.html"), res = 300)
save_as_image(ft2, path = here::here("plots/TabS3.png"), res = 300)
