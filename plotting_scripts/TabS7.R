#'

library(tidyverse)
library(flextable)
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

tb_tidied3 <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_perm.csv"))
tb_tidied5 <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_posthoc.csv"))
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: ")) %>%
        str_remove(fixed(": value")) %>%
        str_replace("exp_id", "strain") %>%
        str_replace("exp_labgroup", "labgroup")
}

ft2 <- tb_tidied3 %>%
    select(ii, gradient, trait_pre, st, term, statistic, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab, ii) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre)
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_at(j = "Model", i = 1:4) %>%
    merge_at(j = "Model", i = 5:8) %>%
    merge_at(j = "Model", i = 9:12) %>%
    merge_at(j = "Model", i = 13:16) %>%
    merge_at(j = "Model", i = 17:20) %>%
    merge_at(j = "Model", i = 21:24) %>%
    merge_v(j = c("Gradient", "Trait")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Trait", "Model", "Term"), align = "center", part = "all") %>%
    #line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", j = c("Term", "Chisq", "P"), i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft2, path = here::here("plots/TabS7.html"), res = 300)
save_as_image(ft2, path = here::here("plots/TabS7.png"), res = 300)
