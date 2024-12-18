#' This script compares the traits in pairs of populations
#' It outputs four tables
#' 1. pairs_anova.csv is the stat table of simple anova
#' 2. pairs_perm.csv is the stat table of permutation test.
#' 3. pairs_perm_obv.csv is the value of simple anova
#' 4. pairs_perm_raw.csv is the raw permutation output used for generating histogram

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
    pivot_longer(cols = -c(1:13), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    left_join(isolates) %>%
    ungroup()

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
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    unnest(mod_tidied) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, statistic, df, p.value) %>%
    mutate(statistic = round(statistic, 2)) %>%
    mutate(ast = map_chr(p.value, turn_p_to_asteriks))

write_csv(tb_tidied, paste0(folder_data, "phenotypes/plants/pairs_anova.csv"))

# 2.2 Permutation ----
tb$chisq_perm <- list(NA)
n_perms = 1000

for (i in 1:nrow(tb)) {
    st <- tb$st[i]
    dat <- tb$dat[[i]]
    list_chisq <- list(NA)
    # Permute data
    for (j in 1:n_perms) {
        set.seed(j)
        d <- dat
        d$population <- d$population[order(runif(nrow(d)))]
        eval(parse(text = st))
        mod_tidied <- Anova(mod, type = 3) %>% tidy()
        list_chisq[[j]] <- select(mod_tidied, term, statistic)
        if (j %% 100 == 0) cat("\n", j)
    }
    tb$chisq_perm[[i]] <- bind_rows(list_chisq)
    cat("\n", i)
}

get_perm_p <- function (chisq_perm, chisq_obv) {
    #' Rank and find the oberved chisq among the permutation values
    temp <- bind_rows(chisq_perm, tibble(statistic = chisq_obv)) %>%
        arrange(desc(statistic)) %>%
        mutate(rank = row_number()) %>%
        filter(statistic == chisq_obv) %>%
        unique()
    return((temp$rank-1)/nrow(chisq_perm))
}
tb_tidied2 <- tb %>%
    left_join(traits) %>%
    mutate(
        ii = 1:n(),
        chisq_obv = map(mod, ~Anova(.x, type = 3) %>% tidy %>% select(term, statistic))
    ) %>%
    unnest(chisq_obv)

tb_tidied2$chisq_perm2 <- NA
for (i in 1:nrow(tb_tidied2)) tb_tidied2$chisq_perm2[i] <- list(filter(tb_tidied2$chisq_perm[[i]], term == tb_tidied2$term[i]))

tb_tidied3 <- tb_tidied2 %>%
    left_join(traits) %>%
    mutate(
        statistic = round(statistic, 2),
        p_value = map2_dbl(chisq_perm2, statistic, get_perm_p),
        siglab = map_chr(p_value, clean_p_lab)
    ) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, statistic, p_value, siglab)

write_csv(tb_tidied3, paste0(folder_data, "phenotypes/plants/pairs_perm.csv"))


# For histogram ----
pairs_perm_raw <- tb_tidied2 %>%
    select(-dat, -mod, -term, -statistic) %>%
    mutate(trait_pre = factor(trait_pre, traits$trait_pre)) %>%
    unnest(chisq_perm) %>%
    filter(!str_detect(term, "Intercept"))
pairs_perm_obv <- tb_tidied2 %>%
    filter(!str_detect(term, "Intercept")) %>%
    select(gradient, trait_pre, trait, term, statistic)
write_csv(pairs_perm_obv, paste0(folder_data, "phenotypes/plants/pairs_perm_obv.csv"))
write_csv(pairs_perm_raw, paste0(folder_data, "phenotypes/plants/pairs_perm_raw.csv"))

