#' This script plots the trait data of plant experiments

library(tidyverse)
library(cowplot)
library(ggh4x)
library(flextable)
library(ggh4x) # for nested facets
library(broom.mixed) # for tidying the model outputs
library(lme4) # for lmer
library(car) # for anova
library(boot) # for bootstrapping
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))
# 1. Prepare data ----
plants_n <- plants %>%
    filter(population != "control", exp_plant == "sativa", gradient == "elevation") %>%
    filter(exp_nitrogen == "without nitrogen") %>%
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

# 3. Run pairwise stat ----
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
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)"
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

# 4. Make the table ----
ft <- tb_tidied %>%
    left_join(traits) %>%
    select(Type = trait_type, Trait = trait_pre, Model = st, Effect = effect, Term = term, Estimate = t0, `95% CIs` = cis, ` ` = signlab, trait_abr) %>%
    # Clean the table
    #filter(Predictor != "(Intercept)") %>%
    mutate(
        Model = str_replace(Model, "value", trait_abr) %>% str_remove("mod <-"),
        Trait = factor(Trait, traits$trait_pre)
    ) %>%
    arrange(Trait) %>%
    flextable() %>%
    valign(valign = "top") %>%
    merge_v(j = 1:4) %>%
    valign(j = 1:3, valign = "center") %>%
    align(j = c("Type", "Trait", "Effect"), align = "center", part = "all") %>%
    hline(i = c(4,7,11,14)) %>%
    autofit() %>%
    width(j = 2, 2) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "population")) %>%
    bg(bg = "grey90", j = "Effect", i = ~Effect == "fixed") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    line_spacing(j = "Model", space = 1.5) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "plants/01-sativa_traits_table.png"), res = 300)
#save_as_html(ft, path = paste0(folder_phenotypes, "plants/03-sativa_nitrogen_table.html"))



# 5. Plot ----
plot_boxes <- function (plants, gra, plant) {
    # gra = "elevation"
    # plant = "sativa"
    strips = strip_nested(
        background_y = elem_list_rect(
            color = NA,
            fill = c("white")
        ),
        text_y = elem_list_text(
            angle = 0
        ),
        bleed = T,
        by_layer_y = F,
        clip = "off", size = "variable"
    )

    plants %>%
        filter(gradient == gra) %>%
        filter(population != "control", exp_plant == plant, exp_nitrogen == "without nitrogen") %>%
        select(-primary_root_nodule_number, -lateral_root_nodule_number) %>%
        select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
        group_by(gradient, population, exp_plant) %>%
        filter(nodule_number <100) %>%
        pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
        left_join(traits) %>%
        arrange(trait_type) %>%
        group_by(gradient, population, exp_plant, trait_type, trait_pre, value) %>%
        count() %>%
        ggplot(aes(x = population, y = value)) +
        #geom_violin(aes(fill = population), position = "identity", alpha = 0.5, color = NA, trim = T) +
        geom_boxplot(aes(fill = population), alpha = 0.5, width = .3, outlier.size = -1) +
        #geom_point(alpha = .8, shape = 16, aes(color = population)) +
        #geom_dotplot(aes(fill = population), binwidth = .1, binaxis = "y", stackgroups = TRUE, binpositions = "all", stackdir = "center") +
        geom_point(alpha = .2, shape = 16, aes(color = population, size = n)) +
        #geom_jitter(size = .5, shape = 21, width = .2) +
        scale_fill_manual(values = population_colors) +
        scale_color_manual(values = population_colors) +
        scale_size_continuous(range = c(1,10)) +
        coord_flip(clip = "off") +
        facet_nested_wrap(trait_type + trait_pre ~., ncol = 1, strip.position = "left", axes = "all", scales = "free", solo_line = T, nest_line = element_line(color = "grey30", linetype = 1, linewidth = 1), strip = strips) +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(color = NA, fill = "grey95"),
            strip.background = element_rect(color = NA),
            legend.position = "top",
            legend.key = element_rect(fill = NA, color = NA),
            legend.key.height = unit(10, "mm"),
            legend.text = element_text(size = 8),
            legend.title = element_blank(),
            legend.background = element_blank(),
            legend.box.margin = unit(c(0,0,-5,0), "mm")
        ) +
        guides(size = "none", fill = guide_legend(override.aes = list(color = NA, size = 0, shape = 0))) +
        labs()
}
p <- plot_boxes(plants, "elevation", "sativa")
ggsave(paste0(folder_phenotypes, "plants/01-sativa_traits.png"), p, width = 8, height = 5)



if (F) {

set.seed(1)
p <- plants %>%
    filter(population != "control", exp_plant == "sativa", exp_nitrogen == "without nitrogen") %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    left_join(traits) %>%
    ggplot(aes(x = population, y = value, fill = population)) +
    geom_boxplot(outlier.shape = -1, alpha = 0.5) +
    geom_jitter(size = .5, shape = 21, width = .2) +
    scale_fill_manual(values = population_colors) +
    coord_cartesian(clip = "off") +
    facet_nested(
        gradient~trait_type+trait_pre, switch = "y", scales = "free", independent = "all", render_empty = F,
        axes = "x", remove_labels = "none",
        strip = strip_nested(bleed=T, clip = "off", size = "variable", text_x = element_text(size = 10), background_x = elem_list_rect(color = NA, fill = c(rep("grey90", 4), rep("white", 10))))) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = "grey90"),
    ) +
    guides(fill = "none") +
    labs(y = "")

}
