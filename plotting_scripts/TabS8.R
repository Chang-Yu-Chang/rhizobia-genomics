#'

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))

pairs_rn_posthoc <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_posthoc.csv"))

clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_replace("exp_id", "strain") %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_remove("glmer\\(|lmer\\(")
}

ft <- pairs_rn_posthoc %>%
    select(ii, gradient, trait_pre, st, temperature, t.ratio, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Temperature = temperature, `T ratio` = t.ratio, P = siglab, ii) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre)
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Trait", "Model", "Temperature", "T ratio"), align = "center", part = "all") %>%
    autofit() %>%
    width(j = "Model", width = 4) %>%
    # Lines and background
    hline(i = c(4,7,11,15,18,22)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "pink", j = "P", i = ~str_detect(P, "\\*")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft, path = here::here("plots/TabS8.html"), res = 300)
save_as_image(ft, path = here::here("plots/TabS8.png"), res = 300)



