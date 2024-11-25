#' This script compares the traits in pairs of populations

library(tidyverse)
library(ggh4x)
library(flextable)
library(broom.mixed) # for tidying the model outputs
library(glmmTMB) # for checking GLMM assumptions
library(DHARMa) # for checking GLMM assumptions
library(lme4) # for lmer
library(car) # for anova
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))


# 1. Prepare data ----
plants_n <- plants %>%
    filter(population != "control", exp_plant == "sativa", exp_nitrogen == "N-") %>%
    select(-nodule_shape, -nodule_size, -nodule_color) %>%
    select(-primary_root_nodule_number, -lateral_root_nodule_number) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:13), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    left_join(isolates) %>%
    ungroup()

# 2. Check assumptions ----
dat <- plants_n %>%
    filter(trait_pre == "shoot height (cm)")
    #filter(trait_pre == "nodule number")
    #filter(trait_pre == "leaf number")
    #filter(trait_pre == "longest petiole\nlength (cm)")

dat %>%
    ggplot() +
    geom_histogram(aes(x = value)) +
    facet_grid(exp_nitrogen~population) +
    theme_bw() +
    theme() +
    guides() +
    labs()

mod <- dat %>%
    filter(gradient == "elevation") %>%
    #filter(exp_nitrogen == "N+") %>%
    #mutate(value = log(value)) %>%
    glmmTMB(value ~ population, data = .)
    #glmmTMB(value ~ population, data = ., family = "nbinom2", ziformula = ~1)
    #glmmTMB(value ~ population * exp_nitrogen, data = ., family = "gaussian")
plot(simulateResiduals(mod))
Anova(mod, type = 3)


# 3. Run model ----
# 3.1 Run ANOVA ----
do_stat <- function (dat, st) {
    d <- dat
    eval(parse(text = st))
    return(mod)
}
tb <- tibble(
    gradient = c(rep("elevation", 4), rep("urbanization", 4)),
    trait_pre = c(
        # Elevation traits
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number",
        # "longest petiole\nlength (cm)",
        # Urbanization traits
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number"
        # "primary root\nlength (cm)",
        # "lateral root number",
        # "longest lateral\nroot length(cm)"
    ),
    st = c(
        # Elevation traits
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        #"mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        # Urbanization traits
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)"
        # "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        # "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        # "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)"
    )
) %>%
    mutate(
        dat = map2(gradient, trait_pre, ~filter(plants_n, gradient ==.x, trait_pre == .y)),
        mod = map2(dat, st, do_stat)
    )

# Tidy up
tb_tidied <- tb %>%
    left_join(traits) %>%
    arrange(gradient, trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    unnest(mod_tidied) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, statistic, df, p.value) %>%
    mutate(statistic = round(statistic, 2)) %>%
    mutate(ast = map(p.value, turn_p_to_asteriks))

# 3.2 Permutation ----
tb$chisq_perm <- list(NA)
n_perms = 1000

for (i in 1:nrow(tb)) {
    st <- tb$st[i]
    dat <- tb$dat[[i]]
    list_chisq <- rep(NA, n_perms)
    # Permute data
    for (j in 1:n_perms) {
        set.seed(j)
        d <- dat
        d$population <- d$population[order(runif(nrow(d)))]
        eval(parse(text = st))
        mod_tidied <- Anova(mod, type = 3) %>% tidy()
        list_chisq[j] <- mod_tidied$statistic[mod_tidied$term == "population"]
        if (j %% 100 == 0) cat("\n", j)
    }
    tb$chisq_perm[[i]] <- tibble(chisq = list_chisq)
    cat("\n", i)
}

get_perm_p <- function (chisq_perm, chisq_obv) {
    #' Rank and find the oberved chisq among the permutation values
    temp <- bind_rows(chisq_perm, tibble(chisq = chisq_obv)) %>%
        arrange(desc(chisq)) %>%
        mutate(rank = row_number()) %>%
        filter(chisq == chisq_obv) %>%
        unique()
    return((temp$rank-1)/nrow(chisq_perm))
}
tb_tidied2 <- tb %>%
    left_join(traits) %>%
    arrange(gradient, trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    mutate(chisq_obv = map_dbl(mod_tidied, ~.x$statistic[.x$term == "population"])) %>%
    # Compute the permutation p value
    mutate(p_value = map2_dbl(chisq_perm, chisq_obv, get_perm_p)) %>%
    relocate(p_value)
    #unnest(c(chisq_obv, chisq_perm)) %>%
    #select(-dat, -mod, -mod_tidied)
    #select(ii, gradient, trait_type, trait)
    #select(ii, gradient, trait_type, trait_pre, st, term, ind, t0, ci_lower, ci_upper) %>%
    # Clean the numeric
    # mutate(across(c(t0, ci_lower, ci_upper), ~round(.x, 2))) %>%
    # mutate(cis = paste0("[", ci_lower, ", ", ci_upper, "]")) %>%
    # mutate(signlab = ifelse(sign(ci_lower) * sign(ci_upper)==T, "*", "n.s."))




# 4. Make table ----
# 4.1 Anova ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: "))
}
ft <- tb_tidied %>%
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

save_as_image(ft, path = paste0(folder_phenotypes, "plants/pairs_tab_lab_anova.png"), res = 300)

# 4.2 Permutation ----
ft2 <- tb_tidied2 %>%
    select(-dat, -mod, -chisq_perm) %>%
    unnest(mod_tidied) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, statistic, df, p_value) %>%
    mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
    #mutate(p_value = map_chr(p_value, edit_p)) %>%
    mutate(siglab = map_chr(p_value, clean_p_lab)) %>%
    mutate(statistic = round(statistic, 2)) %>%
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

save_as_image(ft2, path = paste0(folder_phenotypes, "plants/pairs_tab_lab_perm.png"), res = 300)

# 5. Plot permutation ----
tb_tidied2_obv <- distinct(tb_tidied2, chisq_obv, .keep_all = T) %>%
    mutate(trait_pre = factor(trait_pre, traits$trait_pre))

p <- tb_tidied2 %>%
    select(-dat, -mod) %>%
    mutate(trait_pre = factor(trait_pre, traits$trait_pre)) %>%
    unnest(chisq_perm) %>%
    ggplot() +
    geom_histogram(aes(x = chisq), color = "black", fill = "white") +
    geom_vline(data = tb_tidied2_obv, aes(xintercept = chisq_obv), color = "red", linetype = 2) +
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
ggsave(paste0(folder_phenotypes, "plants/pairs_tab_lab_perm_hist.png"), p, width = 6, height = 8)
