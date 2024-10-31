#' This script plots the reaction norm of nitrogen treatmetns

library(tidyverse)
library(cowplot)
library(flextable)
library(ggh4x)
library(broom.mixed)
library(lme4)
library(car)
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))


# Reaction norm for N supplement
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
    left_join(traits)

plants_n_summ <- plants_n %>%
    group_by(exp_nitrogen, trait_type, trait_pre, population) %>%
    summarize(trait_mean = mean(value), trait_sem = sd(value)/sqrt(n())) %>%
    mutate(lower = trait_mean-qnorm(0.975)*trait_sem, upper = trait_mean+qnorm(0.975)*trait_sem)

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

ft <- tb %>%
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

ft

#save_as_html(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.html"))
save_as_image(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.png"), res = 200)


#
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




