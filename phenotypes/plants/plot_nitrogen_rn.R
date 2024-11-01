#' This script plots the reaction norm of nitrogen treatments
#' 1. Prepare the table
#' 2. Check model assumptions
#' 3. Run models
#' 4. Table
#' 5. Plot

library(tidyverse)
library(cowplot)
library(flextable)
library(ggh4x) # for nested facets
library(broom.mixed) # for tidying the model outputs
library(glmmTMB) # for checking GLMM assumptions
library(DHARMa) # for checking GLMM assumptions
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
    filter(population != "control", exp_plant == "sativa", gradient == "elevation") %>%
    mutate(exp_nitrogen = case_when(
        exp_nitrogen == "without nitrogen" ~ "N-",
        exp_nitrogen == "with nitrogen" ~ "N+"
    )) %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    select(-primary_root_nodule_number, -lateral_root_nodule_number) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    ungroup()

plants_n_summ <- plants_n %>%
    group_by(exp_nitrogen, trait_type, trait_pre, population) %>%
    summarize(trait_mean = mean(value), trait_sem = sd(value)/sqrt(n())) %>%
    mutate(lower = trait_mean-qnorm(0.975)*trait_sem, upper = trait_mean+qnorm(0.975)*trait_sem)


# 2. Check assumptions ----
dat <- plants_n %>%
    #filter(trait_pre == "shoot height (cm)")
    #filter(trait_pre == "nodule number")
    #filter(trait_pre == "leaf number")
    filter(trait_pre == "longest petiole\nlength (cm)")
# xx <- dat %>%
#     filter(exp_nitrogen == "N-", population == "low elevation") %>%
#     pull(value)
# ks.test(xx, function (x) pnorm(x, mean = mean(xx), sd = sd(xx))) # test for gaussian
# ks.test(log(xx), function (x) pnorm(x, mean = mean(log(xx)), sd = sd(log(xx)))) # test for log normal
# ks.test(xx, function (x) ppois(x, mean(xx))) # test for poisson
#ks.test(xx, function (x) rbinom(x, mean(xx))) # test for poisson

dat %>%
    ggplot() +
    geom_histogram(aes(x = value)) +
    facet_grid(exp_nitrogen~population) +
    theme_bw() +
    theme() +
    guides() +
    labs()

mod <- dat %>%
    #filter(exp_nitrogen == "N+") %>%
    #mutate(value = log(value)) %>%
    #glmmTMB(value ~ population, data = ., family = "nbinom2", ziformula = ~1)
    glmmTMB(value ~ population * exp_nitrogen, data = ., family = "gaussian")
plot(simulateResiduals(mod))
Anova(mod, type = 3)


# 3. Fit the model ----
tidy_mod <- function (mod) {
    x <- tidy(mod)$estimate
    return(x)
}
boot_fun <- function(data, indices) {
    #' Define the bootstrapping function. This is generic depending on the exact model `mod` and the tidy function `tidy_mod`
    d <- data[indices, ]; eval(parse(text = st))
    #mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)
    return(tidy_mod(mod)) # Collect the fixed effect estimates
}
get_boot_cis <- function (boot_result) {
    #' Get the CIs of bootstrap value
    funn <- function(x) {
        tb <- as_tibble(x[4][[1]])
        colnames(tb) <- c("ci_lev", "ci_lower", "ci_upper")
        return(tb)
    }
    tibble(
        #term = tidy(mod)$term[1:4],
        t0 = boot_result$t0,
        ind = 1:length(boot_result$t0)
    ) %>%
        drop_na(t0) %>%
        mutate(
            conf = map(ind, ~boot.ci(boot_result, index = .x, type = "norm")),
            ci = map(conf, funn)
        ) %>%
        unnest(ci) %>%
        select(-conf) %>%
        select(ind, t0, everything())
}
do_stat <- function (dat, st) {
    d <- dat
    eval(parse(text = st))
    return(mod)
}

## Bootstrap
set.seed(1)
tb <- tibble(
    trait_pre = rep(c(
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number",
        "longest petiole\nlength (cm)"
    ), each = 1),
    st = rep(c(
        # "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_id) + (0+exp_nitrogen|exp_id), data = d)",
        # "mod <- lmer(value ~ population*exp_nitrogen + (exp_nitrogen|exp_id), data = d)",
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)",

        # "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_id) + (0+exp_nitrogen|exp_id), family = 'poisson', data = d)",
        # "mod <- glmer(value ~ population*exp_nitrogen + (exp_nitrogen|exp_id), family = 'poisson', data = d)",
        "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), family = 'poisson', data = d)",

        # "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_id) + (0+exp_nitrogen|exp_id), data = d)",
        # "mod <- lmer(value ~ population*exp_nitrogen + (exp_nitrogen|exp_id), data = d)",
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)",

        # "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_id) + (0+exp_nitrogen|exp_id), family = 'poisson', data = d)",
        # "mod <- glmer(value ~ population*exp_nitrogen + (exp_nitrogen|exp_id), family = 'poisson', data = d)",
        "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), family = 'poisson', data = d)",

        # "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_id) + (0+exp_nitrogen|exp_id), data = d)",
        # "mod <- lmer(value ~ population*exp_nitrogen + (exp_nitrogen|exp_id), data = d)",
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)"
    ), each = 1)
) %>%
    mutate(
        dat = map(trait_pre, ~filter(plants_n, trait_pre == .x)),
        mod = map2(dat, st, do_stat)
    )

tb$mod_boot <- list(NA)
tb$mod_cis <- list(NA)
for (i in 1:nrow(tb)) {
    st <- tb$st[i]
    dat <- tb$dat[[i]]
    tb$mod_boot[[i]] <- boot(data = dat, statistic = boot_fun, R = 100)
    tb$mod_cis[[i]] <- get_boot_cis(tb$mod_boot[[i]])
}

tb_tidied <- tb %>%
    mutate(mod_tidied = map(mod, tidy)) %>%
    unnest(c(mod_tidied, mod_cis)) %>%
    select(-dat, -mod_boot, mod) %>%
    select(trait_pre, st, ind, t0, ci_lower, ci_upper, effect, group, term, estimate, statistic) %>%
    mutate(across(c(t0, ci_lower, ci_upper), ~round(.x, 2))) %>%
    mutate(cis = paste0("[", ci_lower, ", ", ci_upper, "]")) %>%
    mutate(signlab = ifelse(sign(ci_lower) * sign(ci_upper)==T, "*", "n.s."))

# Simple test ANOVA
tb_tidied2 <- tb %>%
    mutate(ano = map(mod, ~tidy(Anova(.x,type=3)))) %>%
    unnest(ano) %>%
    select(trait_pre, st, term, statistic, df, p.value)


# 4. Make the table ----
ft <- tb_tidied %>%
    select(Trait = trait_pre, Model = st, Effect = effect, Term = term, Estimate = t0, `95% CIs` = cis, ` ` = signlab) %>%
    # Clean the table
    #filter(Predictor != "(Intercept)") %>%
    mutate(
        Model = str_replace(Model, "value", "trait") %>% str_remove("mod <-"),
        Trait = factor(Trait, traits$trait_pre)
    ) %>%
    arrange(Trait) %>%
    flextable() %>%
    valign(valign = "top") %>%
    merge_v(j = 1:3) %>%
    valign(j = 1:3, valign = "center") %>%
    align(j = c("Trait", "Effect"), align = "center", part = "all") %>%
    hline(i = c(6,11,17,22)) %>%
    autofit() %>%
    width(j = 2, 3) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~Term == "population1:exp_nitrogen1") %>%
    bg(bg = "grey90", j = "Effect", i = ~Effect == "fixed") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    line_spacing(j = "Model", space = 1.5) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.png"), res = 300)
#save_as_html(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.html"))

ft2 <- tb_tidied2 %>%
    left_join(traits) %>%
    mutate(statistic = round(statistic, 2)) %>%
    select(Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df =df,  p.value) %>%
    # Clean the table
    #filter(Predictor != "(Intercept)") %>%
    mutate(
        Model = str_replace(Model, "value", "trait") %>% str_remove("mod <-"),
        Trait = factor(Trait, traits$trait_pre),
        P = map_chr(p.value, clean_p_lab)
        #siglab = map_chr(p.value, detect_sig)
    ) %>%
    select(-p.value) %>%
    arrange(Trait) %>%
    flextable() %>%
    autofit() %>%
    merge_v(j = 1:3) %>%
    width(j = "Model", 3) %>%
    valign(j = 1:3, valign = "center") %>%
    align(j = c("Trait", "Term"), align = "center", part = "all") %>%
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()


save_as_image(ft2, path = paste0(folder_phenotypes, "plants/03b-sativa_nitrogen_table.png"), res = 300)














# 5. Plot the reaction norm ----
p <- plants_n %>%
    group_by(gradient, population, exp_plant, exp_nitrogen, trait_type, trait_pre, value) %>%
    count() %>%
    ggplot(aes(x = exp_nitrogen, y = value)) +
    #geom_boxplot(aes(fill = population), outlier.shape = -1, alpha = 0.5) +
    #geom_violin(aes(fill = population), position = "identity", alpha = 0.5, color = NA, trim = T) +
    #geom_dotplot(aes(x = exp_nitrogen, y = value, color = population, fill = population), binaxis = 'y', stackdir = 'center', dotsize = .1) +
    geom_point(aes(color = population, size = n), alpha = .5, shape = 16, position = position_dodge(width = .3)) +
    geom_linerange(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, ymin = lower, ymax = upper),linewidth = 1, position = position_dodge2(width = .3)) +
    geom_line(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, group = population), linewidth = 1, position = position_dodge(width = .3)) +
    geom_point(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean), size = 2, shape = 21, stroke = 1, fill = "white", position = position_dodge2(width = .3)) +
    scale_color_manual(values = population_colors) +
    scale_size_continuous(range = c(1, 15)) +
    facet_nested(
        ~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none",
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("grey90", 3), rep("white", 7))))) +
    theme_bw() +
    theme(
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = "grey90")
    ) +
    guides(
        color = guide_legend(title = "Population", override.aes = list(size = 5)),
        size = guide_legend(title = "Sample size")
    ) +
    labs(x = "Nitrogen treatment", y = "")

ggsave(paste0(folder_phenotypes, "plants/03-sativa_nitrogen.png"), p, width = 12, height = 5)







if (F) {
    # tr <- "shoot height (cm)"
    # st <- "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)"
    # # Fit a GLMM
    # dat <- plants_n %>% filter(trait_pre == tr)
    # mod <- lmer(as.formula(ff), data = dat)
    # tidy_mod(mod)
    # # Bootstrap the model
    # set.seed(1)
    # boot_results <- boot(data = dat, statistic = boot_fun, R = 1000)
    # get_boot_cis(boot_results)

}

if (F) {
    # Old code

    tb <- tibble(
        trait_pre = rep(c(
            "shoot height (cm)",
            "nodule number",
            "leaf color",
            "leaf number",
            "longest petiole\nlength (cm)"
        ), each = 4),
        ff = rep(c(
            "value ~ population * exp_nitrogen + (1|exp_id)",
            "value ~ population * exp_nitrogen + (1|exp_id) + (0 + exp_nitrogen | exp_id)",
            "value ~ population * exp_nitrogen + (exp_nitrogen | exp_id)",
            "value ~ population * exp_nitrogen + (1 | exp_nitrogen : exp_id)"
        ), 5),
        fam = rep(c("gaussian", "poisson", "gaussian", "poisson", "gaussian"), each = 4)
    ) %>%
        mutate(
            dat = map(trait_pre, ~filter(plants_n, trait_pre == .x)),
            mod = pmap(list(x=dat, y=ff, z=fam), function(x,y,z) glmer(as.formula(y), data = x, family = z)),
            mod_tided = map(mod, ~tidy(Anova(.x, type=3)))
        ) %>%
        unnest(cols = mod_tided) %>%
        select(trait_pre, ff, term, fam, statistic, df, p.value) %>%
        filter(ff != "value ~ population * exp_nitrogen + (1|exp_id)")

    ft <- tb_tidied %>%
        select(Trait = trait_pre, Family = fam, Formula = ff, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
        # Clean the table
        filter(Predictor != "(Intercept)") %>%
        mutate(
            Formula = str_replace(Formula, "value", Trait),
            Chisq = round(Chisq, 2),
            `Signif.` = map_chr(P, detect_sig),
            P = map_chr(P, edit_p),
            Trait = factor(Trait, traits$trait_pre)
        ) %>%
        arrange(Trait) %>%
        flextable() %>%
        valign(valign = "top") %>%
        merge_v(j = 1:3) %>%
        valign(j = 1:3, valign = "center") %>%
        hline(i = seq(3, nrow_part(.), 3)) %>%
        autofit() %>%
        bg(bg = "white", part = "all") %>%
        bg(bg = "grey90", i = seq(3, nrow_part(.), 3), j = 4:8) %>%
        style(part = "header", pr_t = fp_text_default(bold = T)) %>%
        style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
        highlight(j = 8, i = ~ `Signif.` != "n.s." & Predictor == "population:exp_nitrogen", color = "yellow") %>%
        fix_border_issues()
}

