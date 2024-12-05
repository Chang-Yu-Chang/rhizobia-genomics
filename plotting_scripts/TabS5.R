#'

library(tidyverse)
library(flextable)
library(broom.mixed) # for tidying the model outputs
library(lme4) # for lmer
library(car) # for anova
library(boot) # for bootstrapping
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

# 1. Prepare the table ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

plants_n <- plants %>%
    filter(nodule_number <100) %>%
    filter(population != "control", exp_plant == "sativa", gradient == "elevation") %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -primary_root_nodule_number, -lateral_root_nodule_number) %>%
    group_by(gradient, population) %>%
    pivot_longer(cols = -c(1:12), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    ungroup()

plants_n_summ <- plants_n %>%
    group_by(exp_nitrogen, trait_type, trait_pre, population) %>%
    summarize(trait_mean = mean(value), trait_sem = sd(value)/sqrt(n())) %>%
    mutate(lower = trait_mean-qnorm(0.975)*trait_sem, upper = trait_mean+qnorm(0.975)*trait_sem)

# 2.1 Anova ----
do_stat <- function (dat, st) {
    d <- dat
    eval(parse(text = st))
    return(mod)
}
tb <- tibble(
    trait_pre = c(
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number"
    ),
    st = c(
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id) + (1|exp_labgroup), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id) + (1|exp_labgroup), data = d)",
        "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id) + (1|exp_labgroup), family = 'poisson', data = d)"
    ),
) %>%
    mutate(
        dat = map(trait_pre, ~filter(plants_n, trait_pre == .x)),
        mod = map2(dat, st, do_stat)
    )

# Tidy up
tb_tidied <- tb %>%
    left_join(traits) %>%
    arrange(trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    unnest(mod_tidied) %>%
    select(ii, trait_type, trait_pre, st, term, statistic, df, p.value) %>%
    mutate(statistic = round(statistic, 2))



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
        list_chisq[j] <- mod_tidied$statistic[mod_tidied$term == "population:exp_nitrogen"]
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
    arrange(trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    mutate(chisq_obv = map_dbl(mod_tidied, ~.x$statistic[.x$term == "population"])) %>%
    # Compute the permutation p value
    mutate(p_value = map2_dbl(chisq_perm, chisq_obv, get_perm_p)) %>%
    relocate(p_value)

# 3. Make the table ----
# 3.1 anova ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: "))
}
ft <- tb_tidied %>%
    select(Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df, p.value, ii) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        P = map_chr(p.value, clean_p_lab),
        P = ifelse(str_detect(Term, "Intercept"), "", P)
    ) %>%
    select(-ii, -p.value) %>%
    arrange(Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Type", "Trait", "Model")) %>%
    valign(j = c("Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

# 4.2 permutation ----
ft2 <- tb_tidied2 %>%
    select(-dat, -mod, -chisq_perm) %>%
    unnest(mod_tidied) %>%
    #select(ii, trait_pre, st, term, statistic, df, p_value) %>%
    mutate(ast = map_chr(p_value, turn_p_to_asteriks)) %>%
    mutate(siglab = map_chr(p_value, clean_p_lab)) %>%
    mutate(statistic = round(statistic, 2)) %>%
    select(ii, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Estimate = statistic, P = siglab) %>%
    # Clean the table
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre)
    ) %>%
    select(-ii) %>%
    arrange(Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Type", "Trait", "Model")) %>%
    valign(j = c("Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()


save_as_html(ft2, path = here::here("plots/TabS5.html"), res = 300)
save_as_image(ft2, path = here::here("plots/TabS5.png"), res = 300)

