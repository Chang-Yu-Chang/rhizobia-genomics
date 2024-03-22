#' This script generates phenotypic comparison for PA populations

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(vegan) # for computing jaccard
source(here::here("analysis/00-metadata.R"))


# Read plant data
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
plants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
plants_long <- read_csv(paste0(folder_data, "temp/23-plants_long.csv"))

# Read growth rate data
gc_prm_summs <- read_csv(paste0(folder_data, 'temp/21-gc_prm_summs.csv'))
isolates_gc <- gc_prm_summs %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    pivot_longer(cols = -c(exp_id, temperature), names_to = "trait") %>%
    unite(trait, trait, temperature) %>%
    left_join(isolates)


# Panel A. Plot he growth trait
compute_trait_mean <- function (isolates_gc, tra = "r_30c", pop = "PA") {
    igcl <- isolates_gc %>%
        filter(trait == tra) %>%
        filter(population == pop)
    igcm <- igcl %>%
        group_by(population, site_group, trait) %>%
        summarize(mean = mean(value, na.rm = T), sem = sd(value, na.rm = T) / sqrt(n()), .groups = "keep")
    return(list(igcl = igcl, igcm = igcm))
}
plot_dots <- function (igcl, igcm) {
    set.seed(1)
    igcl %>%
        ggplot() +
        geom_rect(data = distinct(igcm, site_group), aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = igcm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = igcm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = unique(igcm$trait))
}
test_sign <- function (p) {
    if (p < 0.001) {
        x <- "***"
    } else if (p < 0.01) {
        x <- "**"
    } else if (p < 0.05) {
        x <- "*"
    } else {
        x <- "n.s."
    }
    return(x)
}
plot_pair <- function (tra = "r_30c", pop = "PA", y_axis = expression(r~at~30*degree*C~(1/hr))) {
    t2_1 <- compute_trait_mean(isolates_gc, tra, pop)
    igcl <- t2_1$igcl
    igcm <- t2_1$igcm
    isolates_test <- filter(isolates_gc, trait == tra, population == pop)
    names(isolates_test)[names(isolates_test) == tra] <- "value"
    mod <- lmer(value ~ site_group + (1|site), data = isolates_test)
    mod2_1 <- Anova(mod, type = 3) # no
    sigs <- tibble(population = factor(pop), sig = test_sign(mod2_1[2,3]))
    min_value <- min(igcl$value, na.rm = T)
    max_value <- max(igcl$value, na.rm = T)
    p2 <- igcl %>%
        ggplot() +
        geom_tile(data = igcm, aes(x = site_group, y = mean, fill = site_group), alpha = 0.2, height = Inf, width = 1) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", width = 0.1, height = 0) +
        geom_point(data = igcm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = igcm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.3) +
        # Significance bars
        annotate("segment", x = 1, xend = 2, y = max_value*1.05, yend = max_value*1.05) +
        geom_text(data = sigs, aes(label = sig), x = 1.5, y = max_value*1.05, vjust = -1) +
        scale_fill_manual(values = site_group_colors) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(limits = c(min_value, max_value*1.1)) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(2, "mm"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            panel.grid.major.y = element_line(color = "grey90", linetype = 1, linewidth = .5),
            panel.grid.minor.y = element_line(color = "grey90", linetype = 2, linewidth = .2),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = y_axis)
    return(p2)
}

unique(isolates_gc$trait)
list_va <- rep(list(NA), 12)
list_va[[1]] <- plot_pair(tra = "r_25c", pop = "PA", y_axis = expression(r~at~25*degree*C~(1/hr)))
list_va[[2]] <- plot_pair(tra = "lag_25c", pop = "PA", y_axis = expression(lag~at~25*degree*C~(hr)))
list_va[[3]] <- plot_pair(tra = "maxOD_25c", pop = "PA", y_axis = expression(maxOD~at~25*degree*C))
list_va[[4]] <- plot_pair(tra = "r_30c", pop = "PA", y_axis = expression(r~at~30*degree*C~(1/hr)))
list_va[[5]] <- plot_pair(tra = "lag_30c", pop = "PA", y_axis = expression(lag~at~30*degree*C~(hr)))
list_va[[6]] <- plot_pair(tra = "maxOD_30c", pop = "PA", y_axis = expression(maxOD~at~30*degree*C))
list_va[[7]] <- plot_pair(tra = "r_35c", pop = "PA", y_axis = expression(r~at~35*degree*C~(1/hr)))
list_va[[8]] <- plot_pair(tra = "lag_35c", pop = "PA", y_axis = expression(lag~at~35*degree*C~(hr)))
list_va[[9]] <- plot_pair(tra = "maxOD_35c", pop = "PA", y_axis = expression(maxOD~at~35*degree*C))


# Plot the symbiosis traits comparing the two populations
compute_trait_mean2 <- function (plants_long, tra = "dry_weight_mg", pop = "PA") {
    pl <- plants_long %>%
        filter(trait == tra) %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots2 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(data = distinct(plm, site_group), aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "shoot biomass (mg)")
}
plot_pair2 <- function (tra = "dry_weight_mg", pop = "PA", y_axis = "dry_weight_mg") {
    t3_1 <- compute_trait_mean2(plants_long, tra = tra, pop = pop)
    pl <- t3_1$pl
    plm <- t3_1$plm
    plants_test <- filter(plants, population == pop, exp_id != "control")
    names(plants_test)[names(plants_test) == tra] <- "value"
    plants_test <- drop_na(plants_test, value)
    mod <- lmer(value ~ site_group + (1|site), data = plants_test)
    mod2_1 <- Anova(mod, type = 3) # no
    sigs <- tibble(population = factor(pop), sig = test_sign(mod2_1[2,3]))
    min_value <- min(pl$value, na.rm = T)
    max_value <- max(pl$value, na.rm = T)

    bar_y=55
    p3 <- pl %>%
        ggplot() +
        geom_tile(data = plm, aes(x = site_group, y = mean, fill = site_group), alpha = 0.2, height = Inf, width = 1) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", width = 0.1, height = 0) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.3) +
        # Significance bars
        annotate("segment", x = 1, xend = 2, y = max_value*1.05, yend = max_value*1.05) +
        geom_text(data = sigs, aes(label = sig), x = 1.5, y = max_value*1.05, vjust = -1) +
        scale_fill_manual(values = site_group_colors) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(limits = c(min_value, max_value*1.1)) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(2, "mm"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            panel.grid.major.y = element_line(color = "grey90", linetype = 1, linewidth = .5),
            panel.grid.minor.y = element_line(color = "grey90", linetype = 2, linewidth = .2),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        guides(fill = "none") +
        labs(x = " ", y = y_axis)
    return(p3)
}


list_va[[10]] <- plot_pair2(tra = "dry_weight_mg", pop = "PA", y_axis = "shoot biomass (mg)")
list_va[[11]] <- plot_pair2(tra = "nodule_number", pop = "PA", y_axis = "# of nodules")
list_va[[12]] <- plot_pair2(tra = "root_weight_mg", pop = "PA", y_axis = "root biomass (mg)")

p <- plot_grid(plotlist = list_va[1:12], nrow = 4, align = "vh", axis = "lrbt", labels = c("A", rep("", 8), "B", "", ""))
ggsave(here::here("plots/FigS7.png"), p, width = 10, height = 12)






