#' This script compares the traits in pairs of populations
#' 1. Prepare data
#' 2. Check assumptions
#' 3. Run models
#' 4. Make the stat tables
#' 5. Plot

library(tidyverse)
library(cowplot)
library(ggh4x)
library(flextable)
library(broom.mixed) # for tidying the model outputs
library(lme4) # for lmer
library(car) # for anova
library(boot) # for bootstrapping
source(here::here("metadata.R"))
options(contrasts=c("contr.sum", "contr.poly"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

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


# 1. Prepare data ----
plants_n <- plants %>%
    filter(population != "control", exp_plant == "sativa", exp_nitrogen == "without nitrogen") %>%
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
    left_join(isolates) %>%
    ungroup()

# 2. Check assumptions ----

# 3. Run models ----
# 3.1 Single anova ----
tb <- tibble(
    gradient = c(rep("elevation", 5), rep("urbanization", 7)),
    trait_pre = c(
        # Elevation traits
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number",
        "longest petiole\nlength (cm)",
        # Urbanization traits
        "shoot height (cm)",
        "nodule number",
        "leaf color",
        "leaf number",
        "primary root\nlength (cm)",
        "lateral root number",
        "longest lateral\nroot length(cm)"
    ),
    st = c(
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)",
        "mod <- glmer(value ~ population + (1|exp_id), family = 'poisson', data = d)",
        "mod <- lmer(value ~ population + (1|exp_id), data = d)"
    )
) %>%
    mutate(
        dat = map2(gradient, trait_pre, ~filter(plants_n, gradient ==.x, trait_pre == .y)),
        mod = map2(dat, st, do_stat)
    )

# Tidy up
tb_tidied <- tb %>%
    left_join(traits) %>%
    arrange(gradient, trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, ~Anova(.x, type = 3) %>% tidy())) %>%
    unnest(mod_tidied) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, statistic, df, p.value) %>%
    mutate(statistic = round(statistic, 2))


# 3.2 Bootstrap ----
set.seed(1)


tb$mod_boot <- list(NA)
tb$mod_cis <- list(NA)
for (i in 1:nrow(tb)) {
    st <- tb$st[i]
    dat <- tb$dat[[i]]
    tb$mod_boot[[i]] <- boot(data = dat, statistic = boot_fun, R = 100)
    tb$mod_cis[[i]] <- get_boot_cis(tb$mod_boot[[i]])
    cat(i)
}

tb_tidied2 <- tb %>%
    left_join(traits) %>%
    arrange(gradient, trait_type, trait_pre) %>%
    mutate(ii = 1:n()) %>%
    mutate(mod_tidied = map(mod, tidy)) %>%
    unnest(c(mod_tidied, mod_cis)) %>%
    select(-dat, -mod_boot, mod) %>%
    left_join(traits) %>%
    select(ii, gradient, trait_type, trait_pre, st, term, ind, t0, ci_lower, ci_upper) %>%
    # Clean the numberic
    mutate(across(c(t0, ci_lower, ci_upper), ~round(.x, 2))) %>%
    mutate(cis = paste0("[", ci_lower, ", ", ci_upper, "]")) %>%
    mutate(signlab = ifelse(sign(ci_lower) * sign(ci_upper)==T, "*", "n.s."))


# 4. Make the tables ----
# 4.1 anova ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: "))
}
ft <- tb_tidied %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df, p.value, ii) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
        P = map_chr(p.value, clean_p_lab),
        P = ifelse(str_detect(Term, "Intercept"), "", P)
    ) %>%
    select(-ii, -p.value) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Type", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "plants/pairs_tab_anova.png"), res = 300)

# 4.2 Bootstrap ----
ft2 <- tb_tidied2 %>%
    select(Gradient = gradient, Type = trait_type, Trait = trait_pre, Model = st, Term = term, Estimate = t0, `95% CIs` = cis, ` ` = signlab, ii) %>%
    # Clean the table
    filter(!str_detect(Term, "sd__")) %>%
    mutate(
        Model = map2_chr(Model, ii, ~clean_model_string(.x,.y)),
        #Model = str_replace(Model, "value", trait_abr) %>% str_remove("mod <-"),
        Trait = factor(Trait, traits$trait_pre),
        ` ` = ifelse(str_detect(Term, "Intercept|sd__"), "", ` `)
    ) %>%
    select(-ii) %>%
    arrange(Gradient, Trait) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Type", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Type", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Type", "Trait", "Term"), align = "center", part = "all") %>%
    line_spacing(j = "Trait", space = 1.5) %>%
    # Lines and background
    hline(i = seq(2, nrow_part(.), 2)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "lightpink", i = ~str_detect(Term, "pop")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft2, path = paste0(folder_phenotypes, "plants/pairs_tab_boot.png"), res = 300)



# 5. Plot ----
# 5.1 Boxplot ----
plot_boxes <- function (plants_n, gra, plant) {
    # gra = "elevation"
    # plant = "sativa"
    strips = strip_nested(
        background_y = elem_list_rect(
            color = NA,
            fill = c("white")
        ),
        text_y = elem_list_text(
            size = 8,
            angle = 0
        ),
        bleed = T,
        by_layer_y = F,
        clip = "off", size = "variable"
    )

    plants_n %>%
        filter(gradient == gra) %>%
        left_join(traits) %>%
        arrange(trait_type) %>%
        group_by(gradient, population, exp_plant, trait_type, trait_pre, value, exp_id) %>%
        count() %>%
        ggplot() +
        geom_boxplot(aes(x = population, y = value, fill = population), alpha = 0.3, width = .7, outlier.size = -1) +
        geom_point(aes(x = population, group = exp_id, y = value, color = population, size = n), alpha = .4, shape = 16, position = position_dodge(width = .7)) +
        scale_fill_manual(values = population_colors) +
        scale_color_manual(values = population_colors) +
        scale_size_continuous(range = c(.5,3)) +
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
p1 <- plot_boxes(plants_n, "elevation", "sativa")
p2 <- plot_boxes(plants_n, "urbanization", "sativa")
p <- plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2], scale = .95) +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    draw_text(c("Elevation", "Urbanization"), x = c(0.05, 0.55), y = 0.98, hjust = 0)
ggsave(paste0(folder_phenotypes, "plants/pairs_boxes.png"), p, width = 8, height = 6)

