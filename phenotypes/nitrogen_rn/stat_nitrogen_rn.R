#'

library(tidyverse)
library(broom.mixed) # for tidying the model outputs
library(lme4) # for lmer
library(car) # for anova
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

# 1. Prepare the table ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

plants_n <- plants %>%
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
    mutate(
        ii = 1:n(),
        chisq_obv = map(mod, ~Anova(.x, type = 3) %>% tidy %>% select(term, statistic))
    ) %>%
    unnest(chisq_obv)

# Decouple the three terms. pop, nit, and pop:nit
tb_tidied2$chisq_perm2 <- NA
for (i in 1:nrow(tb_tidied2)) tb_tidied2$chisq_perm2[i] <- list(filter(tb_tidied2$chisq_perm[[i]], term == tb_tidied2$term[i]))

tb_tidied3 <- tb_tidied2 %>%
    left_join(traits) %>%
    mutate(
        statistic = round(statistic, 2),
        p_value = map2_dbl(chisq_perm2, statistic, get_perm_p),
        siglab = map_chr(p_value, clean_p_lab)
    ) %>%
    select(ii, trait_type, trait_pre, st, term, statistic, p_value, siglab)

write_csv(tb_tidied3, paste0(folder_data, "phenotypes/nitrogen_rn/nitrogen_rn_perm.csv"))
