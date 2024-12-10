#' This script compares the growth traits in pairs of populations

library(tidyverse)
library(ggh4x)
library(flextable)
library(emmeans) # for computing estimated marginal means
library(broom.mixed) # for tidying the model outputs
# library(glmmTMB) # for checking GLMM assumptions
# library(DHARMa) # for checking GLMM assumptions
library(lme4) # for lmer
library(car) # for anova
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

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
do_emm <- function (mod) {
    emmeans(mod, specs = "population", by = "temperature") %>% pairs() %>% as_tibble()
}
conditonal_filter <- function (x, y){
    if (y %in% "lag") {
        filter(gtwl, gradient == x, trait == y) %>% filter(temperature != "40c")
    } else if (y %in% c("r", "maxOD")) filter(gtwl, gradient == x, trait == y)
}

tb <- tibble(
    gradient = rep(c("elevation", "urbanization"), each = 3),
    trait = rep(c("r", "lag", "maxOD"), 2),
    st = rep(c(
        "mod <- lmer(value ~ population*temperature + (1|temperature:exp_id), data = d)",
        "mod <- lmer(value ~ population*temperature + (1|temperature:exp_id), data = d)",
        "mod <- lmer(value ~ population*temperature + (1|temperature:exp_id), data = d)"
    ), 2)
) %>%
    mutate(
        dat = map2(gradient, trait, conditonal_filter),
        mod = map2(dat, st, do_stat)
    )

# Tidy up
tb_tidied <- tb %>%
    left_join(traits) %>%
    arrange(gradient, trait_pre) %>%
    mutate(
        ii = 1:n(),
        mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy()),
        mod_emm = map(mod, do_emm)
    ) %>%
    select(-trait_abr, -trait_pre2)



# 3.2 Permutation ----
tb$chisq_perm <- list(NA)
tb$emm_perm <- list(NA)
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
        list_emm[[j]] <- do_emm(mod) %>% select(temperature, t.ratio)
        if (j %% 100 == 0) cat("\n", j)
    }
    tb$chisq_perm[[i]] <- bind_rows(list_chisq)
    tb$emm_perm[[i]] <- bind_rows(list_emm)
    cat("\n", i)
}

get_perm_p <- function (chisq_perm, chisq_obv) {
    #' Rank and find the oberved chisq among the permutation values
    temp <- bind_rows(chisq_perm, tibble(statistic = chisq_obv)) %>%
        arrange(desc(statistic)) %>%
        #group_by(term) %>%
        mutate(rank = row_number()) %>%
        #left_join(rename(chisq_obv, chisq = chisq_obv)) %>%
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

# Decouple the three terms. Pop, temp, and pop:temp
tb_tidied2$chisq_perm2 <- NA
for (i in 1:nrow(tb_tidied2)) tb_tidied2$chisq_perm2[i] <- list(filter(tb_tidied2$chisq_perm[[i]], term == tb_tidied2$term[i]))

tb_tidied3 <- tb_tidied2 %>%
    left_join(traits) %>%
    mutate(
        statistic = round(statistic, 2),
        p_value = map2_dbl(chisq_perm2, statistic, get_perm_p),
        siglab = map_chr(p_value, clean_p_lab)
    ) %>%
    select(ii, gradient, trait_pre, st, term, statistic, p_value, siglab)

# 3.3 Post hoc tuckey comparison ----
tb_tidied4 <- tb %>%
    mutate(
        ii = 1:n(),
        emm_obv = map(mod, ~select(do_emm(.x), temperature, t.ratio))
    ) %>%
    unnest(emm_obv)

get_perm_p2 <- function (emm_perm, emm_obv) {
    #' Rank and find the oberved chisq among the permutation values
    temp <- bind_rows(emm_perm, tibble(t.ratio = emm_obv)) %>%
        arrange(desc(t.ratio)) %>%
        mutate(rank = row_number()) %>%
        filter(t.ratio == emm_obv) %>%
        unique()
    p_value <- (temp$rank-1)/nrow(emm_perm)
    if (p_value > 0.5) p_value <- 1-p_value
    return(p_value)
}

tb_tidied5 <- tb_tidied4 %>%
    left_join(traits) %>%
    mutate(
        t.ratio = round(t.ratio, 2),
        p_value = map2_dbl(emm_perm, t.ratio, get_perm_p2),
        siglab = map_chr(p_value, clean_p_lab)
    ) %>%
    select(ii, gradient, trait_pre, st, temperature, t.ratio, p_value, siglab)

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
    unnest(mod_tidied) %>%
    mutate(
        statistic = round(statistic, 2),
        ast = map_chr(p.value, turn_p_to_asteriks),
        st = map2_chr(st, ii, ~clean_model_string(.x,.y)),
        trait = factor(trait, traits$trait),
        p.value = map_chr(p.value, clean_p_lab)
    ) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df, P = p.value) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "growth/pairs_rn_anova.png"), res = 300)

# 4.2 Permutation ----
ft2 <- tb_tidied3 %>%
    select(ii, gradient, trait_pre, st, term, statistic, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P_perm = siglab, ii) %>%
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
    valign(j = c("Gradient", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft2, path = paste0(folder_phenotypes, "growth/pairs_rn_perm.png"), res = 300)
write_csv(tb_tidied3, paste0(folder_phenotypes, "growth/pairs_rn_perm.csv"))

# 4.3 Post hoc tukey test ----
ft3 <- tb_tidied5 %>%
    select(ii, gradient, trait_pre, st, temperature, t.ratio, p_value, siglab) %>%
    select(Gradient = gradient, Trait = trait_pre, Model = st, Temperature = temperature, T_ratio = t.ratio, P_perm = siglab, ii) %>%
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
    valign(j = c("Gradient", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Trait", "Temperature"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = ~str_detect(Temperature, "40c")) %>%
    hline(i = ~str_detect(Temperature, "35c") & str_detect(Trait, "lag")) %>%
    bg(bg = "white", part = "all") %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft3, path = paste0(folder_phenotypes, "growth/pairs_rn_posthoc.png"), res = 300)
write_csv(tb_tidied5, paste0(folder_phenotypes, "growth/pairs_rn_posthoc.csv"))






