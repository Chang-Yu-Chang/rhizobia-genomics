#' This script generates phenotypic comparison for all traits

renv::load()
library(tidyverse)
library(cowplot)
library(janitor)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
#library(vegan) # for computing jaccard
source(here::here("metadata.R"))

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


isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv")) %>%
    clean_names() %>%
    filter(temperature != "40c") %>%
    select(exp_id, r, lag, max_od, temperature) %>%
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, max_od)) %>%
    drop_na %>%
    left_join(isolates)

# Planel A. plot the growth trait
set.seed(1)
get_p <- function (tra, pop) {
    tb <- gts %>%
        mutate(population = factor(population, c("VA", "PA"))) %>%
        filter(population == pop) %>%
        rename(dep = {{tra}})


    mod <- lmer(dep ~ site_group + (1|site), data = tb)
    p_value <- Anova(mod, type = 3)[2,3]
    vmax <- max(`[`(tb, "dep"))*.95
    return(tibble(vmax=vmax, p_value=p_value))
}


tb_stat <- tibble(trait = colnames(gts)[2:10]) %>%
    rowwise() %>%
    mutate(mod = list(get_p(trait, "VA"))) %>%
    unnest(mod) %>%
    rowwise() %>%
    mutate(sig = detect_sig(p_value)) %>%
    mutate(trait = str_replace(trait, "max_od", "maxod")) %>%
    separate(col = trait, into = c("trait_group", "temperature")) %>%
    mutate(trait_group = factor(trait_group, c("r", "lag", "maxod"))) %>%
    group_by(trait_group) %>%
    mutate(vmax = max(vmax))

tb <- gts %>%
    filter(population == "VA") %>%
    pivot_longer(ends_with("c"), names_to = "trait") %>%
    mutate(trait = str_replace(trait, "max_od", "maxod")) %>%
    separate(col = trait, into = c("trait_group", "temperature")) %>%
    mutate(trait_group = factor(trait_group, c("r", "lag", "maxod")))

plot_box <- function (tra = "r", yy) {
    tb %>%
        filter(trait_group == tra) %>%
        ggplot() +
        geom_boxplot(aes(x = site_group, y = value, fill = site_group)) +
        geom_jitter(aes(x = site_group, y = value, color = site_group), size = 2, width = .1, shape = 21, stroke = 1) +
        # Stat
        geom_segment(data = filter(tb_stat, trait_group == tra), x = 1.05, xend = 1.95, aes(y = vmax, yend = vmax)) +
        geom_text(data = filter(tb_stat, trait_group == tra), x = 1.5, aes(y = vmax*1.05, label = sig)) +
        scale_fill_manual(values = alpha(site_group_colors, 0.5)) +
        scale_color_manual(values = alpha(site_group_colors, 1)) +
        facet_grid(~ temperature, scales = "free_y") +
        theme_bw() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 10),
            panel.border = element_rect(color = "grey10", fill = NA),
            strip.background = element_blank(),
            plot.background = element_blank(),
            legend.title = element_blank(),
            legend.position = "top",
            legend.box.margin = margin(0,0,0,0),
            legend.margin = margin(0,0,0,0),

        ) +
        guides() +
        labs(y = yy)
}
p1 <- plot_box(tra = "r", yy = "growth rate (1/hr)")
p2 <- plot_box(tra = "lag", yy = "lag time (hr)") + theme(legend.position = "none")
p3 <- plot_box(tra = "maxod", yy = "yield") + theme(legend.position = "none")

# Plot the growth trait


p <- plot_grid(p1,p2,p3, ncol = 1, align = "hv", axis = "blr", rel_heights = c(1.2,1,1)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS3.png"), p, width = 6, height = 8)

if (F) {

# Panel A. Plot he growth trait
compute_trait_mean <- function (isolates_gc, tra = "r_30c", pop = "VA") {
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
plot_pair <- function (tra = "r_30c", pop = "VA", y_axis = expression(r~at~30*degree*C~(1/hr))) {
    t2_1 <- compute_trait_mean(gts, tra, pop)
    igcl <- t2_1$igcl
    igcm <- t2_1$igcm
    isolates_test <- filter(gts, trait == tra, population == pop)
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

unique(gts$trait)
list_va <- rep(list(NA), 12)
list_va[[1]] <- plot_pair(tra = "r_25c", pop = "VA", y_axis = expression(r[T==25](1/hr)))
list_va[[2]] <- plot_pair(tra = "lag_25c", pop = "VA", y_axis = expression(t[T==25](hr)))
list_va[[3]] <- plot_pair(tra = "maxOD_25c", pop = "VA", y_axis = expression(x[T==25]))
list_va[[4]] <- plot_pair(tra = "r_30c", pop = "VA", y_axis = expression(r[T==30](1/hr)))
list_va[[5]] <- plot_pair(tra = "lag_30c", pop = "VA", y_axis = expression(t[T==30](hr)))
list_va[[6]] <- plot_pair(tra = "maxOD_30c", pop = "VA", y_axis = expression(x[T==30]))
list_va[[7]] <- plot_pair(tra = "r_35c", pop = "VA", y_axis = expression(r[T==35](1/hr)))
list_va[[8]] <- plot_pair(tra = "lag_35c", pop = "VA", y_axis = expression(t[T==35](hr)))
list_va[[9]] <- plot_pair(tra = "maxOD_35c", pop = "VA", y_axis = expression(x[T==35]))


# Plot the symbiosis traits comparing the two populations
compute_trait_mean2 <- function (plants_long, tra = "dry_weight_mg", pop = "VA") {
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
plot_pair2 <- function (tra = "r_30c", pop = "VA", y_axis = "dry_weight_mg") {
    t3_1 <- compute_trait_mean2(plants_long,tra = tra, pop = pop)
    pl <- t3_1$pl
    plm <- t3_1$plm
    plants_test <- filter(plants, population == pop, exp_id != "control")
    names(plants_test)[names(plants_test) == tra] <- "value"
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

list_va[[10]] <- plot_pair2(tra = "shoot_biomass_mg", pop = "VA", y_axis = "shoot biomass (mg)")
list_va[[11]] <- plot_pair2(tra = "nodule_count", pop = "VA", y_axis = "# of nodules")
list_va[[12]] <- plot_pair2(tra = "root_biomass_mg", pop = "VA", y_axis = "root biomass (mg)")

p <- plot_grid(plotlist = list_va[1:12], nrow = 4, align = "vh", axis = "lrbt", labels = c("A", rep("", 8), "B", "", ""))
ggsave(here::here("plots/FigS3.png"), p, width = 10, height = 12)







}

