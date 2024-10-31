#' This script plots the reaction norm of nitrogen treatmetns

library(tidyverse)
library(janitor)
library(cowplot)
library(flextable)
library(ggh4x) # for nested facets
library(broom.mixed) # for cleanup the model output
library(lme4) # for LMM
library(car) # For anova
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))


isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    rename(yield = maxOD) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

# reaction norm
tb <- tibble(
    gradient = rep(c("elevation", "urbanization"), each = 3*4),
    trait = rep(rep(c("r", "lag", "yield"), each = 4), 2),
    ff = rep(c(
        "value ~ population * temperature + (1|exp_id) + (1|well)",
        "value ~ population * temperature + (1|exp_id) + (0 + temperature | exp_id) + (1|well)",
        "value ~ population * temperature + (temperature | exp_id) + (1|well)",
        "value ~ population * temperature + (1 | temperature : exp_id) + (1|well)"
        ), 6)
) %>%
    #filter(!(temperature == "40c" & trait == "lag")) %>%
    mutate(
        dat = map2(gradient, trait, ~filter(gtwl, gradient == .x, trait == .y)),
        mod = map2(dat, ff, ~lmer(as.formula(.y), data = .x)),
        mod_tided = map(mod, ~tidy(Anova(.x, type=3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    select(gradient, trait, ff, term, statistic, df, p.value) %>%
    filter(ff != "value ~ population * temperature + (1|exp_id) + (1|well)")


ft <- tb %>%
    select(Gradient = gradient, Trait = trait, Formula = ff, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
    # Clean the table
    filter(Predictor != "(Intercept)") %>%
    mutate(
        Formula = str_replace(Formula, "value", Trait),
        Chisq = round(Chisq, 2),
        `Signif.` = map_chr(P, detect_sig),
        P = map_chr(P, edit_p),
        Trait = factor(Trait, c("r", "lag", "yield"))
    ) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    valign(valign = "top") %>%
    merge_v(j = 1:3) %>%
    valign(j = 1:3, valign = "center") %>%
    hline(i = seq(3, nrow_part(.), 3)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = seq(3, nrow_part(.), 3), j = 4:8) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    highlight(j = 8, i = ~ `Signif.` != "n.s." & Predictor == "population:temperature", color = "yellow") %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_data, "phenotypes/growth/03-rn_table.png"), res = 200)



plants

if (F) {

# Compute the trait ~ population X temperature
tb_poptemp <- tibble(
    grad = rep(c("elevation", "urbanization"), each = 3),
    res = rep(c("r", "lag", "yield"), 2),
    ff = rep(
        c("r ~ population*temperature + (1|site) + (1|genome_id)",
          "lag ~ population*temperature + (1|site) + (1|genome_id)",
          "max_od ~ population*temperature + (1|site) + (1|genome_id)"
        ), 2
    )
) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gtw, gradient == grad)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    rowwise() %>%
    mutate(signif = detect_sig(p.value))




# Reaction norm
gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(temperature, well, exp_id, r, lag, max_od) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    mutate(trait = case_when(
        trait == "r" ~ "growth rate (1/hr)",
        trait == "lag" ~ "lag time (hr)",
        trait == "max_od" ~ "yield [OD]"
    ))

tb_pertemp <- tibble(
    grad = rep(c("elevation", "urbanization"), each = 12),
    res = rep(rep(c("r", "lag", "yield"), each = 4), 2),
    temp = rep(c("25c", "30c", "35c", "40c"), 6),
    ff = rep(rep(
        c("r ~ population + (1|site)",
          "lag ~ population + (1|site)",
          "max_od ~ population + (1|site)"
        ), each = 4), 2
    )
) %>%
    filter(!(temp == "40c" & res == "lag")) %>%
    rowwise() %>%
    mutate(
        dat = list(filter(gtw, gradient == grad, temperature == temp)),
        mod = list(lmer(as.formula(ff), data = dat)),
        mod_tided = list(tidy(Anova(mod, type = 3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    rowwise() %>%
    select(trait_pre, ff, term, fam, statistic, df, p.value)
    mutate(signif = detect_sig(p.value))

# Compute the mean
gtwlm <- gtwl %>%
    group_by(gradient, temperature, trait, population) %>%
    summarize(mean_value = mean(value, na.rm = T), ci_value = qnorm(0.975) * sd(value, na.rm = T) / sqrt(sum(!is.na(value))), n = sum(!is.na(value))) %>%
    group_by(gradient, temperature, trait) %>%
    mutate(max_mean_value = max(mean_value, na.rm = T))
tb_pertemp <- tb_pertemp %>%
    filter(term == "population") %>%
    mutate(trait = case_when(
        res == "r" ~ "growth rate (1/hr)",
        res == "lag" ~ "lag time (hr)",
        res == "yield" ~ "yield [OD]"
    )) %>%
    left_join(distinct(select(gtwlm, grad = gradient, trait, temp = temperature, max_mean_value))) %>%
    ungroup()

tb_poptemp <- tb_poptemp %>%
    filter(term == "population:temperature") %>%
    mutate(trait = case_when(
        res == "r" ~ "growth rate (1/hr)",
        res == "lag" ~ "lag time (hr)",
        res == "yield" ~ "yield [OD]"
    ), edited_p = edit_p(p.value)
    ) %>%
    ungroup()
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
        geom_text(data = filter(tb_pertemp, grad == gradf), aes(x = temp, y = max_mean_value, label = signif), vjust = -3, size = 3) +
        # stats temp X population
        geom_text(data = filter(tb_poptemp, grad == gradf), aes(label = paste0("P(pop:temp): ", edited_p)), x = Inf, y = Inf, vjust = 1.1, hjust = 1, size = 3) +
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
            legend.position = "none",
            plot.background = element_blank()
        ) +
        guides() +
        labs(x = expression(Temperature*degree*C))
}
p_rn1 <- plot_rn(gtwl, "elevation")
p_rn2 <- plot_rn(gtwl, "urbanization")

p <- plot_grid(p_rn1, p_rn2, nrow = 1) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_data, "phenotypes/growth/03-rn_withp.png"), p, width = 6, height = 6)
}








# stat
tb <- tibble(
    trait_pre = c("shoot height (cm)", "nodule number", "leaf color", "leaf number", "longest petiole\nlength (cm)"),
    ff = rep("value ~ population*exp_nitrogen + (1|exp_id)", 5),
    fam = c("gaussian", "poisson", "gaussian", "poisson", "gaussian")
) %>%
    mutate(
        dat = map(trait_pre, ~filter(plants_n, trait_pre == .x)),
        mod = pmap(list(x=dat, y=ff, z=fam), function(x,y,z) glmer(as.formula(y), data = x, family = z)),
        mod_tided = map(mod, ~tidy(Anova(.x, type=3)))
    ) %>%
    unnest(cols = mod_tided) %>%
    select(trait_pre, ff, term, fam, statistic, df, p.value)

ft <- tb %>%
    select(Trait = trait_pre, Formula = ff, Family = fam, Predictor = term, Chisq = statistic, df = df, P = p.value) %>%
    rowwise() %>%
    # Clean the table
    mutate(Formula = str_replace(Formula, "value", Trait)) %>%
    filter(Predictor != "(Intercept)") %>%
    mutate(Chisq = round(Chisq, 2), P = edit_p(P)) %>%
    mutate(`Signif.` = detect_sig(P)) %>%
    mutate(Trait = factor(Trait, traits$trait_pre)) %>%
    arrange(Trait) %>%
    #
    flextable() %>%
    valign(valign = "top") %>%
    merge_v(j = 1:3) %>%
    valign(j = 1:3, valign = "center") %>%
    hline(i = c(3,6,9,12)) %>%
    autofit() %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = seq(1,flextable::nrow_part(.),2), j = 4:8) %>%
    #align(j = 1:8, align = "center", part = "all") %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    style(part = "header", j = c("df", "P"), pr_t = fp_text_default(italic = T, bold = T)) %>%
    fix_border_issues()

#save_as_html(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.html"))
save_as_image(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.png"), res = 200)

# plants_n %>%
#     #filter(trait_pre == "shoot height (cm)") %>%
#     filter(trait_pre == "nodule number") %>%
#     #filter(trait_pre == "longest petiole\nlength (cm)") %>%
#     glmer(value ~ exp_nitrogen*population + (1|exp_id), data = ., family = "gaussian") %>%
#     Anova(type = 3)
# # Fit model ("poisson" or "nbinom1")
# modelname = glmmTMB(response~predictor + (1|random), data=data, family="poisson")
# # check assumptions
# plot(simulateResiduals(modelname))
# # if fails, refit model with different error dist
# # if passes, test significance
# # Test significance
# Anova(modelname, type=3)

