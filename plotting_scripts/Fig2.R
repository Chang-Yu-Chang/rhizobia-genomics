#' This script

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(broom.mixed)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for permanova
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv")) %>% slice(1:32)
gcs <- read_csv(paste0(folder_data, "phenotypes/growth/gcs.csv"))
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv"))

set.seed(1)
edit_p <- function (pv) {
    if (pv < 0.001) {
        return("<0.001")
    } else {
        return(as.character(round(pv, 3)))
    }
}
detect_sig <- function (pv) {
    if (pv > 0.05) {
        return("n.s.")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv < 0.001) {
        return("***")
    }
}

# Correct the naming ----
isolates <- isolates %>%
    mutate(gradient = case_when(
        population == "VA" ~ "elevation",
        population == "PA" ~ "urbanization"
    ), .keep = "unused") %>%
    rename(population = site_group)

gtwl <- gtw %>%
    replace_na(list(maxOD = 0)) %>%
    mutate(temperature = factor(temperature, c("25c", "30c", "35c", "40c"))) %>%
    select(-t.r, -startOD) %>%
    pivot_longer(-c(temperature, well, exp_id), names_to = "trait") %>%
    left_join(distinct(isolates, exp_id, .keep_all = T))

# PCA at 30C ----
gts1 <- gtwl %>%
    filter(gradient == "elevation", temperature == "30c") %>%
    pivot_wider(names_from = trait, values_from = value) %>%
    select(population, genome_id, r, lag, maxOD)
gts2 <- gtwl %>%
    filter(gradient == "urbanization", temperature == "30c") %>%
    pivot_wider(names_from = trait, values_from = value) %>%
    select(population, genome_id, r, lag, maxOD)

pca_results <- list(elevation = prcomp(gts1[,-c(1,2)], scale. = TRUE),
                    urbanization = prcomp(gts2[,-c(1,2)], scale. = TRUE))
get_pcvar <- function (pca_result) summary(pca_result)$importance[2, ] %>% round(3) * 100

pcs <- bind_rows(
    as_tibble(pca_results$elevation$x) %>% mutate(population = gts1$population, genome_id = gts1$genome_id),
    as_tibble(pca_results$urbanization$x) %>% mutate(population = gts2$population, genome_id = gts2$genome_id)
) %>%
    left_join(distinct(isolates, gradient, population)) %>%
    mutate(gradient = factor(gradient, c("elevation", "urbanization")))

plot_pca <- function (pcs, grad) {
    pcsi <- pcs %>% filter(gradient == grad)
    dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
    mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = genome_id))
    pcsi %>%
        ggplot() +
        geom_point(aes(x = PC1, y = PC2, color = population), shape = 21, stroke = 1, size = 2) +
        stat_ellipse(aes(x = PC1, y = PC2, fill = population), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 3))) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        scale_color_manual(values = site_group_colors) +
        scale_fill_manual(values = site_group_colors, name = "population") +
        scale_x_continuous(breaks = seq(-8,8,2)) +
        scale_y_continuous(breaks = seq(-8,8,2)) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            plot.background = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            aspect.ratio = 1
        ) +
        guides(color = "none") +
        labs(x = paste0("PC1 (", get_pcvar(pca_results[grad][[1]])[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_results[grad][[1]])[2], "%)"))

}
p_pca1 <- plot_pca(pcs, "elevation")
p_pca2 <- plot_pca(pcs, "urbanization")

p_pca <- plot_grid(p_pca1, p_pca2, nrow = 1, axis = "lr", align = "v", labels = c("B", "C"))

# Stats
pcsi <- pcs %>% filter(gradient == "elevation")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = genome_id))

pcsi <- pcs %>% filter(gradient == "urbanization")
dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = genome_id))

# Compute the per-temp trait ~ population ----
gtw <- read_csv(paste0(folder_data, "phenotypes/growth/gtw.csv")) %>%
    clean_names() %>%
    replace_na(list(max_od = 0)) %>%
    select(exp_id, r, lag, max_od, temperature, well) %>%
    left_join(isolates)

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
    mutate(signif = detect_sig(p.value))


# Compute the trait ~ population X temperature ----
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




# Reaction norm ----
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
        scale_color_manual(values = site_group_colors, name = "population") +
        scale_fill_manual(values = site_group_colors, name = "population") +
        scale_x_discrete(breaks = c("25c", "30c", "35c", "40c"), labels = c(25, 30, 35, 40)) +
        facet_wrap(~trait, scales = "free_y", nrow = 1, strip.position = "left") +
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
#p_rn <- plot_grid(p_rn1, NULL, p_rn2, nrow = 1, axis = "lr", align = "v", rel_widths = c(1,.1,1), labels = c("D", "", "E"))


# Combine figures ----

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig2.png"), scale = 1) +
    #draw_plot(p_pca, x = .3, y = .48, width = .65, height = .43) +
    draw_plot(p_pca1, x = .3, y = .48 , width = .32, height = .43) +
    draw_plot(p_pca2, x = .7, y = .48, width = .3, height = .41) +
    #draw_plot(p_rn, x = .05, y = .05, width = .9, height = .4) +
    draw_plot(p_rn1, x = .03, y = .05, width = .44, height = .4) +
    draw_plot(p_rn2, x = .55, y = .05, width = .43, height = .4) +
    #draw_plot(p_sativa, x = .35, y = -0.01, width = .6, height = .32) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 6)
