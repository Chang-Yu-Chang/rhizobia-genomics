#' This script compares the growth traits in pairs of populations

library(tidyverse)
library(emmeans) # for computing estimated marginal means
library(broom.mixed) # for tidying the model outputs
library(lme4) # for lmer
library(car) # for anova
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))
set.seed(1)

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

# 1. Prepare data ----
gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    mutate(trait = factor(trait, c("r", "lag", "maxOD"))) %>%
    left_join(distinct(isolates, exp_id, .keep_all = T))


# 3. Run model ----
# 3.1 Run ANOVA ----
do_stat <- function (dat, st) {
    d <- dat
    eval(parse(text = st))
    return(mod)
}
conditonal_filter <- function (x, y, z){
    if (y %in% "lag") {
        filter(gtwl, gradient == x, trait == y, temperature == z) %>% filter(temperature != "40c")
    } else if (y %in% c("r", "maxOD")) filter(gtwl, gradient == x, trait == y, temperature == z)
}

tb <- tibble(
    gradient = rep(c("elevation", "urbanization"), each = 12),
    temperature = rep(c("25c", "30c", "35c", "40c"), 6),
    trait = rep(rep(c("r", "lag", "maxOD"), each = 4), 2),
    st = rep("mod <- lmer(value ~ population + (1|exp_id), data = d)", 24)
) %>%
    filter(!(temperature == "40c" & trait == "lag")) %>%
    mutate(
        dat = pmap(list(gradient, trait, temperature), conditonal_filter),
        mod = map2(dat, st, do_stat)
    )

# Tidy up
tb_tidied <- tb %>%
    left_join(traits) %>%
    arrange(gradient, trait_pre) %>%
    mutate(
        ii = 1:n(),
        mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())
    ) %>%
    select(ii, gradient, temperature, trait_pre, st, mod_tidied) %>%
    unnest(mod_tidied) %>%
    mutate(
        statistic = round(statistic, 2),
        p.value = round(p.value, 3),
        siglab = map_chr(p.value, clean_p_lab)
    )

write_csv(tb_tidied, paste0(folder_phenotypes, "growth/pairs_anova.csv"))

# 3.2 Permutation ----
tb$chisq_perm <- list(NA)

n_perms = 1000
for (i in 1:nrow(tb)) {
    st <- tb$st[i]
    dat <- tb$dat[[i]]
    list_chisq <- list(NA)
    list_emm <- list(NA)
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
    select(ii, gradient, temperature, trait_pre, st, term, statistic, p_value, siglab)

write_csv(tb_tidied3, paste0(folder_phenotypes, "growth/pairs_perm.csv"))
