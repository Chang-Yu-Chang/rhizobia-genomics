#'

library(tidyverse)
library(flextable)
library(ggh4x)


pairs_lab_anova <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_lab_anova.csv"))
pairs_lab_perm <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_lab_perm.csv"))

pairs_lab_perm_obv <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_lab_perm_obv.csv"))
pairs_lab_perm_forhist <- read_csv(paste0(folder_data, "phenotypes/plants/pairs_lab_perm_forhist.csv"))


# 4.1 Anova ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: "))
}
ft <- pairs_lab_anova %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df, p.value, ii) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        P = map_chr(p.value, clean_p_lab),
        P = ifelse(str_detect(Term, "Intercept"), "", P)
    ) %>%
    select(-ii, -p.value) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Type", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "plants/pairs_lab_anova.png"), res = 300)

# 4.2 Permutation ----
ft2 <- pairs_lab_perm %>%
    # Remove the intercept estimates
    mutate(p_value = ifelse(str_detect(term, "Intercept"), "", p_value)) %>%
    mutate(statistic = ifelse(str_detect(term, "Intercept"), "", statistic)) %>%
    #mutate(siglab = paste0(p_value, " ",ast)) %>%
    #
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Estimate = statistic, P = siglab, ii) %>%
    # Clean the table
    #filter(!str_detect(Term, "sd__")) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        #Model = str_replace(Model, "value", trait_abr) %>% str_remove("mod <-"),
        Trait = factor(Trait, traits$trait_pre),
        P = ifelse(str_detect(Term, "Intercept|sd__"), "", P)
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Type", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft2, path = paste0(folder_phenotypes, "plants/pairs_lab_perm.png"), res = 300)

# 5. Plot permutation ----

p <- pairs_lab_perm_forhist %>%
    ggplot() +
    geom_histogram(aes(x = statistic), color = "black", fill = "white") +
    geom_vline(data = pairs_lab_perm_obv, aes(xintercept = statistic), color = "red", linetype = 2) +
    scale_y_continuous(position = "right") +
    facet_grid2(trait_pre ~ gradient, scale = "free_x", switch = "y", strip = strip_themed(text_y = element_text(angle = 0))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = NA, fill = NA),
        strip.placement = "outside",
        strip.clip = "off"
    ) +
    guides() +
    labs()

ggsave(paste0(folder_phenotypes, "plants/pairs_lab_perm_hist.png"), p, width = 6, height = 6)
