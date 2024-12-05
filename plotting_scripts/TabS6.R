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

# Post hoc tukey test
ft3 <- tb_tidied5 %>%
    select(ii, gradient, trait_pre, st, temperature, t.ratio, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Temperature = temperature, T_ratio = t.ratio, P = siglab, ii) %>%
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
    merge_at(j = "Model", i = 5:7) %>%
    merge_at(j = "Model", i = 8:11) %>%
    merge_at(j = "Model", i = 12:15) %>%
    merge_at(j = "Model", i = 16:18) %>%
    merge_at(j = "Model", i = 19:22) %>%
    merge_v(j = c("Gradient", "Trait")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "top") %>%
    align(j = c("Gradient", "Trait", "Model", "Temperature", "T_ratio"), align = "center", part = "all") %>%
    autofit() %>%
    width(j = "Model", width = 4) %>%
    line_spacing(j = "Model", space = 1.5) %>%
    # Lines and background
    hline(i = c(4,7,11,15,18,22)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "pink", j = "P", i = ~str_detect(P, "\\*")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft3, path = here::here("plots/TabS6.html"), res = 300)
save_as_image(ft3, path = here::here("plots/TabS6.png"), res = 300)



