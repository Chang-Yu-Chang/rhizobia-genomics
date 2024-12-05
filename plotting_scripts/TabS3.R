#' This script compares the traits in pairs of populations inoculated to M. sativa

library(tidyverse)
library(ggh4x)
library(flextable)
library(broom.mixed) # for tidying the model outputs
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

# 2. Run models ----
# 2.1 Run ANOVA ----
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
        # Urbanization traits
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number"
    ),
    st = c(
        # Elevation traits
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        # Urbanization traits
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id) + (1|exp_labgroup), family = 'poisson', data = d)"
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

# 2.2 Permutation ----
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

# 3. Make table ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: ")) %>%
        str_remove(fixed(": value")) %>%
        str_replace("exp_id", "strain") %>%
        str_replace("exp_labgroup", "labgroup")
}
ft2 <- tb_tidied2 %>%
    select(-dat, -mod, -chisq_perm) %>%
    unnest(mod_tidied) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, statistic, df, p_value) %>%
    mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
    mutate(siglab = map_chr(p_value, clean_p_lab)) %>%
    mutate(statistic = round(statistic, 2)) %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Estimate = statistic, P = siglab, ii) %>%
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
    align(j = c("Gradient", "Type", "Trait", "Model", "Term", "P"), align = "center", part = "all") %>%
    autofit() %>%
    width(j = "Model", width = 4) %>%
    line_spacing(j = "Model", space = 1.5) %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_html(ft2, path = here::here("plots/TabS3.html"), res = 300)
save_as_image(ft2, path = here::here("plots/TabS3.png"), res = 300)
