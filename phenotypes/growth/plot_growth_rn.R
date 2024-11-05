#' This script plots the reaction norm of thermal adaptation
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
#options(contrasts=c("contr.sum", "contr.poly"))

# 1. Prepare the data ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    #rename(yield = maxOD) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    drop_na(value)

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

# 2. Check assumption ----
dat <- gtwl %>%
    filter(gradient == "elevation") %>%
    filter(trait == "r")

dat %>%
    ggplot() +
    geom_histogram(aes(x = value)) +
    facet_grid(temperature~population) +
    theme_bw() +
    theme() +
    guides() +
    labs()

mod <- dat %>%
    #filter(exp_nitrogen == "N+") %>%
    #mutate(value = log(value)) %>%
    #glmmTMB(value ~ population, data = ., family = "nbinom2", ziformula = ~1)
    glmmTMB(value ~ population * temperature, data = ., family = "gaussian")
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

conditonal_filter <- function (x, y){
    if (y %in% "lag") {
        filter(gtwl, gradient == x, trait == y) %>% filter(temperature != "40c")
    } else if (y %in% c("r", "maxOD")) filter(gtwl, gradient == x, trait == y)
}

# gtwl <- gtwl %>%
#     filter(temperature != "40c")
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


# gtwl %>%
#     filter(gradient == "elevation", trait == "lag") %>%
#     view
#     lmer(value ~ population*temperature + (1|temperature:exp_id), data = .)
#
# tb$mod[[1]] %>% Anova(type=3)
#
# mod <- tb$mod[[2]]
# tidy(mod)$estimate
# tb$mod[[3]]

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
    select(gradient, trait, st, ind, t0, ci_lower, ci_upper, effect, group, term, estimate, statistic) %>%
    mutate(across(c(t0, ci_lower, ci_upper), ~round(.x, 2))) %>%
    mutate(cis = paste0("[", ci_lower, ", ", ci_upper, "]")) %>%
    mutate(signlab = ifelse(sign(ci_lower) * sign(ci_upper)==T, "*", "n.s."))


# 4. Make the table
ft <- tb_tidied %>%
    select(Gradient = gradient, Trait = trait, Model = st, Effect = effect, Term = term, Estimate = t0, `95% CIs` = cis, ` ` = signlab) %>%
    # Clean the table
    filter(Term != "(Intercept)") %>%
    mutate(
        Model = str_replace(Model, "value", Trait) %>% str_remove("mod <- "),
        Trait = factor(Trait, c("r", "lag", "maxOD"))
    ) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    valign(valign = "top") %>%
    merge_v(j = 1:3) %>%
    valign(j = 1:3, valign = "center") %>%
    #hline(i = seq(9, nrow_part(.), 9)) %>%
    hline(i = c(9,16,25,34,41)) %>%
    width(j = "Model", 2) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "\\:temperature")) %>%
    #bg(bg = "grey90", i = seq(6, nrow_part(.), 6), j = 4:8) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    #highlight(j = 8, i = ~ `Signif.` != "n.s." & Predictor == "population:temperature", color = "yellow") %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_data, "phenotypes/growth/03-rn_table.png"), res = 200)

# 5. Plot reaction norm ----
gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, well, exp_id, r, lag, maxOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "maxOD" ~ "yield [OD]"
    ))

# Compute the mean
gtwlm <- gtwl %>%
    group_by(gradient, temperature, trait, population) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value))) %>%
    group_by(gradient, temperature, trait) %>%
    mutate(max_mean_value = max(mean_value, na.rm = T))
plot_rn <- function (gg, gradf) {
    tt <- filter(gtwlm, gradient == gradf)

    gg %>%
        filter(gradient == gradf) %>%
        ggplot() +
        # Individual replicates
        geom_line(aes(x = temperature, y = value, group = well, color = population), alpha = 0.1) +
        # Mean value
        geom_ribbon(data = tt, aes(x = temperature, ymin = mean_value-ci_value, ymax =  mean_value+ci_value, fill = population, group = population), inherit.aes = FALSE, alpha = 0.2) +
        geom_point(data = tt, aes(x = temperature, y = mean_value, color = population, group = population)) +
        geom_line(data = tt, aes(x = temperature, y = mean_value, color = population, group = population)) +
        # stats per temperature
        #geom_text(data = filter(tb_pertemp, grad == gradf), aes(x = temp, y = max_mean_value, label = signif), vjust = -3, size = 3) +
        # stats temp X population
        #geom_text(data = filter(tb_poptemp, grad == gradf), aes(label = paste0("P(pop:temp): ", edited_p)), x = Inf, y = Inf, vjust = 1.1, hjust = 1, size = 3) +
        scale_color_manual(values = population_colors, name = "population") +
        scale_fill_manual(values = population_colors, name = "population") +
        scale_x_discrete(breaks = c("25c", "30c", "35c", "40c"), labels = c(25, 30, 35, 40)) +
        facet_wrap(~trait, scales = "free_y", ncol = 1, strip.position = "left") +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 10),
            strip.text.y = element_text(size = 10),
            strip.placement = "outside",
            axis.title.x = element_text(size = 10),
            axis.title.y = element_blank(),
            legend.position = "inside",
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            plot.background = element_blank()
        ) +
        guides(fill = guide_legend(override.aes = list(color = NA))) +
        labs(x = expression(paste("Temperature (", degree, "C)")))
}
p_rn1 <- plot_rn(gtwl, "elevation") + theme(legend.position.inside = c(.8, .95))
p_rn2 <- plot_rn(gtwl, "urbanization") + theme(legend.position.inside = c(.2, .95))

p <- plot_grid(
    p_rn1, p_rn2,
    nrow = 1, labels = c("A", "B"), scale = .95
) +
    theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "phenotypes/growth/03-rn_withp.png"), p, width = 8, height = 8)
