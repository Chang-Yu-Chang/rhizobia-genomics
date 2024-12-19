#' This script plots the reaction norm of thermal adaptation

library(tidyverse)
library(cowplot)
library(ggh4x)
library(flextable)
source(here::here("metadata.R"))

pairs_rn_anova <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_anova.csv"))
pairs_rn_perm <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_perm.csv"))
pairs_rn_posthoc <- read_csv(paste0(folder_phenotypes, "growth/pairs_rn_posthoc.csv"))

# Anova table ----
clean_model_string <- function (mod_st, ii) {
    mod_st %>%
        str_replace("value", paste0("trait", as.character(ii))) %>%
        str_remove(fixed("mod <- ")) %>%
        str_remove(fixed(", data = d)")) %>%
        str_replace(fixed("lmer("), fixed("lmer: "))
}
ft <- pairs_rn_anova %>%
    select(ii, Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, df, P = siglab) %>%
    mutate(
        Model = map2_chr(st, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
    ) %>%
    arrange(Gradient, Trait) %>%
    select(-ii) %>%
    flextable() %>%
    autofit() %>%
    # Align and spacing
    merge_v(j = c("Gradient", "Trait", "Model")) %>%
    valign(j = c("Gradient", "Trait", "Model"), valign = "center") %>%
    align(j = c("Gradient", "Trait", "Term"), align = "center", part = "all") %>%
    # Lines and background
    hline(i = seq(4, nrow_part(.), 4)) %>%
    bg(bg = "white", part = "all") %>%
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft, path = paste0(folder_phenotypes, "growth/pairs_rn_anova.png"), res = 300)

# Permutation table ----
ft2 <- pairs_rn_perm %>%
    select(ii, Gradient = gradient, Trait = trait_pre, Model = st, Term = term, Chisq = statistic, P = siglab) %>%
    mutate(
        Model = map2_chr(st, ii, ~clean_model_string(.x,.y)),
        Trait = factor(Trait, traits$trait_pre),
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
    bg(bg = "grey90", i = ~str_detect(Term, ":")) %>%
    style(part = "header", pr_t = fp_text_default(bold = T)) %>%
    fix_border_issues()

save_as_image(ft2, path = paste0(folder_phenotypes, "growth/pairs_rn_perm.png"), res = 300)

# Post hoc tukey test ----
ft3 <- pairs_rn_posthoc %>%
    #select(ii, gradient, trait_pre, st, temperature, t.ratio, p_value, siglab) %>%
    select(ii, Gradient = gradient, Trait = trait_pre, Model = st, Temperature = temperature, `T ratio` = t.ratio, P_perm = siglab, ii) %>%
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

# Plot reaction norm ----
isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

isolates <- isolates %>%
    arrange(population) %>%
    mutate(genome_id = factor(genome_id))

gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T)) %>%
    drop_na(value)

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
ggsave(paste0(folder_data, "phenotypes/growth/pairs_rn.png"), p, width = 8, height = 8)
