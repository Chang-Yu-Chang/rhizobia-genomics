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

# 1. Prepare the table ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

plants_n <- plants %>%
    filter(population != "control", exp_plant == "sativa", gradient == "elevation") %>%
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
# dat <- plants_n %>%
#     #filter(trait_pre == "shoot height (cm)")
#     #filter(trait_pre == "nodule number")
#     #filter(trait_pre == "leaf number")
#     filter(trait_pre == "longest petiole\nlength (cm)")
# xx <- dat %>%
#     filter(exp_nitrogen == "N-", population == "low elevation") %>%
#     pull(value)
# ks.test(xx, function (x) pnorm(x, mean = mean(xx), sd = sd(xx))) # test for gaussian
# ks.test(log(xx), function (x) pnorm(x, mean = mean(log(xx)), sd = sd(log(xx)))) # test for log normal
# ks.test(xx, function (x) ppois(x, mean(xx))) # test for poisson
#ks.test(xx, function (x) rbinom(x, mean(xx))) # test for poisson
#
# dat %>%
#     ggplot() +
#     geom_histogram(aes(x = value)) +
#     facet_grid(exp_nitrogen~population) +
#     theme_bw() +
#     theme() +
#     guides() +
#     labs()
#
# mod <- dat %>%
#     #filter(exp_nitrogen == "N+") %>%
#     #mutate(value = log(value)) %>%
#     #glmmTMB(value ~ population, data = ., family = "nbinom2", ziformula = ~1)
#     glmmTMB(value ~ population * exp_nitrogen, data = ., family = "gaussian")
# plot(simulateResiduals(mod))
# Anova(mod, type = 3)


# 3. Fit the model ----
# 3.1 Anova ----
tb <- tibble(
    trait_pre = rep(c(
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number",
        "longest petiole\nlength (cm)"
    ), each = 1),
    st = rep(c(
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)",
        "mod <- glmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population*exp_nitrogen + (1|exp_nitrogen:exp_id), data = d)",
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

# Tidy up
tb_tidied <- tb %>%
    left_join(traits) %>%
    arrange(trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    unnest(mod_tidied) %>%
    select(ii, trait_type, trait_pre, st, term, statistic, df, p.value) %>%
    mutate(statistic = round(statistic, 2))


# 3.2 Bootstrap ----
## Bootstrap
set.seed(1)

tb$mod_boot <- list(NA)
tb$mod_cis <- list(NA)
for (i in 1:nrow(tb)) {
    st <- tb$st[i]
    dat <- tb$dat[[i]]
    tb$mod_boot[[i]] <- boot(data = dat, statistic = boot_fun, R = 100)
    tb$mod_cis[[i]] <- get_boot_cis(tb$mod_boot[[i]])
}

tb_tidied2 <- tb %>%
    left_join(traits) %>%
    arrange(trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, tidy)) %>%
    unnest(c(mod_tidied, mod_cis)) %>%
    select(-dat, -mod_boot, mod) %>%
    left_join(traits) %>%
    select(ii, trait_type, trait_pre, st, term, ind, t0, ci_lower, ci_upper) %>%
    # Clean the numberic
    mutate(across(c(t0, ci_lower, ci_upper), ~round(.x, 2))) %>%
    mutate(cis = paste0("[", ci_lower, ", ", ci_upper, "]")) %>%
    mutate(signlab = ifelse(sign(ci_lower) * sign(ci_upper)==T, "*", "n.s."))



# 4. Make the table ----
# 4.1 anova ----
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
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "nitrogen_rn/rn_tab_anova.png"), res = 300)


# 4.2 Bootstrap ----
ft2 <- tb_tidied2 %>%
    select(Type = trait_type, Trait = trait_pre, Model = st, Term = term, Estimate = t0, `95% CIs` = cis, ` ` = signlab, ii) %>%
    # Clean the table
    filter(!str_detect(Term, "sd__")) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        ` ` = ifelse(str_detect(Term, "Intercept|sd__"), "", ` `)
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
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

#save_as_image(ft2, path = paste0(folder_phenotypes, "nitrogen_rn/rn_tab_boot.png"), res = 300)

# 5. Plot the reaction norm ----
p <- plants_n %>%
    group_by(gradient, population, exp_plant, exp_nitrogen, trait_type, trait_pre, value) %>%
    count() %>%
    ggplot(aes(x = exp_nitrogen, y = value)) +
    geom_point(aes(color = population, size = n), alpha = .5, shape = 16, position = position_dodge(width = .3)) +
    geom_linerange(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, ymin = lower, ymax = upper), linewidth = 1, position = position_dodge2(width = .3)) +
    geom_line(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean, group = population), linewidth = 1, position = position_dodge(width = .3)) +
    geom_point(data = plants_n_summ, aes(x = exp_nitrogen, y = trait_mean), size = 2, shape = 21, stroke = 1, fill = "white", position = position_dodge2(width = .3)) +
    scale_color_manual(values = population_colors) +
    scale_size_continuous(range = c(.5, 10)) +
    facet_nested(
        ~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none", nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), solo_line = T,
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("white", 3), rep("white", 7))))) +
    theme_bw() +
    theme(
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = NA)
    ) +
    guides(
        color = guide_legend(title = "population", override.aes = list(size = 5)),
        size = guide_legend(title = "sample size")
    ) +
    labs(x = "Nitrogen treatment", y = "")

ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 4)


