#' This script analyses the figures shown in the paper

library(tidyverse)
library(cowplot)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(factoextra) # for plotting pca eclipse
library(flextable) # for making table
library(Hmisc) # for pairwise correlation test
source(here::here("analysis/00-metadata.R"))

# growth curves
gc <- read_csv(paste0(folder_data, 'temp/04-gc.csv'), show_col_types = F)
gc_summ <- read_csv(paste0(folder_data, 'temp/04-gc_summ.csv'), show_col_types = F)
gc.prm <- read_csv(paste0(folder_data, 'temp/04-gc_prm.csv'), show_col_types = F)
gc.prm.stat <- read_csv(paste0(folder_data, 'temp/04-gc_prm_summ.csv'), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/02-isolates_RDP.csv"), show_col_types = F) %>%
    rename(strain = ExpID) %>%
    filter(Genus == "Ensifer", str_sub(strain, 1,1) %in% c("H","L")) %>%
    mutate(strain_site = str_sub(strain, 1, 2), strain_site_group = str_sub(strain, 1, 1))

# mutualism
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/11-treatments_long.csv"), show_col_types = F)
treatments_scaled <- read_csv(paste0(folder_data, "temp/11-treatments_scaled.csv"), show_col_types = F)
treatments_scaled_long <- read_csv(paste0(folder_data, "temp/11-treatments_scaled_long.csv"), show_col_types = F)
# Clean up data
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))


# 0. combine the data from treatments and gc ----
# Subset only Ensifer -----
subset_ensifer <- function(tb) {
    tb %>%
        left_join(select(isolates_RDP, strain, Genus)) %>%
        drop_na()
}

gc <- gc %>% subset_ensifer()
gc_summ <- gc_summ %>% subset_ensifer()
gc.prm <- gc.prm %>% subset_ensifer()
gc.prm.stat <- gc.prm.stat %>% subset_ensifer()

# Combine the data from rhizobia ----
treatments_long_stat <- treatments_long %>%
    filter(trait %in% c("dry_weight_mg", "nodule_number", "root_weight_mg", "total_root_length_px", "branching_frequency_per_px", "network_area_px2", "average_diameter_px")) %>%
    group_by(strain_site_group, strain, trait) %>%
    drop_na(strain, value) %>%
    summarize(trait_mean = mean(value, na.rm = T), trait_sd = sd(value, na.rm = T), n = n())

t1 <- gc.prm.stat %>%
    filter(strain %in% rhizobia_strains) %>%
    select(strain, r, lag, maxOD) %>%
    pivot_longer(cols = -strain, names_to = "trait", values_to = "trait_mean")
t2 <- gc.prm.stat %>%
    filter(strain %in% rhizobia_strains) %>%
    select(strain, r=r.sem, lag=lag.sem, maxOD=maxOD.sem) %>%
    pivot_longer(cols = -strain, names_to = "trait", values_to = "trait_sd")

gc_long_stat <- left_join(t1, t2) %>% mutate(n = 4) %>% mutate(strain_site_group = str_sub(strain, 1, 1))

trait_long_stat <- bind_rows(mutate(gc_long_stat, trait_type = "growth"), mutate(treatments_long_stat, trait_type = "mutualism")) %>%
    select(strain_site_group, strain, everything()) %>%
    arrange(strain_site_group, strain, trait)


# 1. plot growth trait mean by site ----
plot_boxplot_pair <- function (tb, ytrait, ylab = "") {
    tb %>%
        ggplot() +
        geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
        geom_boxplot(aes(x = strain_site_group, y = {{ytrait}}), fill = "white", outlier.size = -1, color = "black") +
        geom_point(aes(x = strain_site_group, y = {{ytrait}}, group = strain, color = strain), shape = 21, size = 2, stroke = 1, fill = NA,
                   position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
        scale_color_manual(values = rep("black", 100)) +
        scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
        #scale_x_discrete(label = c("high elevation", "low elevation")) +
        facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high\nelevation", L="low\nelevation"))) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            #panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
            strip.background = element_rect(color = NA, fill = NA),
            strip.text = element_text(size = 10, color = "black"),
            axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_blank(),
            #axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
            legend.position = "none"
        ) +
        guides(color = "none") +
        labs(x = "", y = ylab)

}
p1 <- plot_boxplot_pair(gc.prm, lag, "lag time (hr)")
p2 <- plot_boxplot_pair(gc.prm, r, expression(growth~rate(h^-1))) #+ theme(legend.position = "top", legend.title = element_blank())
p3 <- plot_boxplot_pair(gc.prm, maxOD, expression(paste("max", "[", OD[600], "]")))

p <- plot_grid(p1, p2, p3, nrow = 1, axis = "tbrl", align = "hv")
ggsave(paste0(folder_data, "temp/22-01-gc_by_site.png"), p, width = 6, height = 4)

## Does the rhizobia sites have effect on any growth trait?
mod <- lmer(lag ~ strain_site_group + (1|strain) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group does not have effect on lag time
mod <- lmer(r ~ strain_site_group + (1|strain) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group does not have effect on r
mod <- lmer(maxOD ~ strain_site_group + (1|strain) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group does not have effect on maxOD

## Does the rhizobia strain have effect on any growth trait?
mod <- lmer(lag ~ strain + (1|strain_site_group) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group has effect on lag time
mod <- lmer(r ~ strain + (1|strain_site_group) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group has effect on r
mod <- lmer(maxOD ~ strain + (1|strain_site_group) + (1|strain_site), data = gc.prm)
Anova(mod, type = 3) # Site group has effect on maxOD

# 1a. plot traits by strains, faceted by sites ----
gc_labels <- gc.prm.stat %>%
    mutate(strain_label = factor(1:n())) %>%
    select(strain, strain_label)
p <- gc.prm %>%
    left_join(gc_labels) %>%
    ggplot() +
    geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
    #geom_boxplot(aes(x = strain, y = lag, fill = strain_site_group), alpha = .6, outlier.size = -1, color = "black") +
    geom_point(aes(x = strain_label, y = lag, group = strain_label, color = strain_label), shape = 21, size = 2, stroke = 1, fill = NA,
               position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
    scale_color_manual(values = rep("black", 100)) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
    #scale_x_discrete(label = c("high", "low")) +
    facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high elevation", L="low elevation"))) +
    theme_classic() +
    theme(
        panel.grid.major.x = element_line(color = "grey80"),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "none"
    ) +
    guides(color = "none") +
    labs(x = "rhizobia strain", y = "lag time (hr)")
#p <- plot_grid(p1, p2, p3, nrow = 1, axis = "tbrl", align = "hv")
ggsave(paste0(folder_data, "temp/22-01a-gc_by_strain.png"), p, width = 4, height = 3)


# 1b. pca of the three growth traits ----
tt <- gc.prm %>%
    select(well, strain_site_group, all_of(c("r", "lag", "maxOD"))) %>%
    drop_na()

pcobj <- tt %>%
    select(-well, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)

p <- fviz_pca_ind(
    pcobj,
    label = "none",
    habillage = tt$strain_site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 16, L = 17), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.2, 0.1),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()
ggsave(paste0(folder_data, "temp/22-01b-pca.png"), p, width = 4, height = 4)



# 2. one extended phenotype trait (root weight) by site group----
p <- treatments %>%
    drop_na(strain, dry_weight_mg) %>%
    plot_boxplot_pair(dry_weight_mg, "above-ground dry weight (mg)") +
    theme(
        axis.text = element_text(size = 10, color = "black")
    )
ggsave(paste0(folder_data, "temp/22-02-weight_by_site.png"), p, width = 2, height = 4)

# 2a. plot one trait by strain, faceted by sites ----
gc_labels <- gc.prm.stat %>%
    mutate(strain_label = factor(1:n())) %>%
    select(strain, strain_label)
p <- treatments_M %>%
    left_join(gc_labels) %>%
    drop_na(strain, dry_weight_mg) %>%
    filter(strain != "control") %>%
    ggplot() +
    geom_rect(data = tibble(strain_site_group = c("H", "L")), aes(fill = strain_site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
    geom_boxplot(aes(x = strain_label, y = dry_weight_mg), fill = "white", outlier.size = -1, color = "black") +
    geom_point(aes(x = strain_label, y = dry_weight_mg, group = strain_label, color = strain_label), shape = 21, size = 2, stroke = 1, fill = NA,
               position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
    scale_color_manual(values = rep("black", 100)) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("high", "low"), breaks = c("H", "L")) +
    facet_grid(~strain_site_group, scales = "free_x", space = "free_x", labeller = labeller(.cols = c(H="high\nelevation", L="low\nelevation"))) +
    theme_classic() +
    theme(
        panel.grid.major.x = element_line(color = "grey80"),
        #panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        panel.spacing.x = unit(0, "mm"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "none"
    ) +
    guides(color = "none") +
    labs(x = "rhizobia strain", y = "above-ground dry weight (mg)")

ggsave(paste0(folder_data, "temp/22-02a-weight_by_strain.png"), p, width = 3, height = 4)


## Does rhizobia strain have effect on dry weight?
mod <- lmer(dry_weight_mg ~ strain + (1|strain_site_group) + (1|plant) + (1|waterblock), data = treatments)
Anova(mod, type = 3) # Site group does not have effect on dry weight


## Does rhizobia sites have effect on dry weight?
mod <- lmer(dry_weight_mg ~ strain_site_group + (1|strain) + (1|plant) + (1|waterblock), data = treatments)
Anova(mod, type = 3) # Site group does not have effect on dry weight


# 3. pca of all extended phenotype traits ----
tt <- treatments_M %>%
    filter(strain != "control") %>%
    select(id, strain_site_group, all_of(traits), -nodule_weight_mg) %>%
    drop_na()

pcobj <- tt %>%
    select(-id, -strain_site_group) %>%
    prcomp(center = TRUE, scale. = TRUE)

p <- fviz_pca_ind(
    pcobj,
    label = "none",
    habillage = tt$strain_site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    scale_shape_manual(values = c(H = 16, L = 17), labels = c("high elevation", "low elevation"), breaks = c("H", "L"), name = "elevation") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.8, 0.9),
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()
ggsave(paste0(folder_data, "temp/22-03-pca.png"), p, width = 4, height = 4)

# 4. pca variables
p <- fviz_pca_var(
    pcobj, repel = T
) +
    scale_x_continuous(limits = c(-1.5,1.5)) +
    scale_y_continuous(limits = c(-1.5,1.5)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    labs()
ggsave(paste0(folder_data, "temp/22-04-pca_var.png"), p, width = 6, height = 6)
# p <- plot_grid(p1, p2, axis = "tblr", align = "h") +
#     theme(plot.background = element_rect(fill = "white", color = NA))

# 5. correlation matrix ----
library(Hmisc) # for pairwise correlation test
calculate_cor <- function (treatments) {
    temp <- treatments %>%
        select(all_of(traits)) %>%
        drop_na(all_of(traits)) %>%
        as.matrix() %>%
        rcorr(type = "spearman")

    tb_cor <- temp$r %>%
        as_tibble() %>%
        mutate(row = colnames(.)) %>%
        pivot_longer(-row, names_to = "col", values_to = "correlation")
    tb_p <- temp$P %>%
        as_tibble() %>%
        mutate(row = colnames(.)) %>%
        pivot_longer(-row, names_to = "col", values_to = "p_value")
    tb <- tb_cor %>%
        left_join(tb_p)
    return(tb)
}
plot_cor_matrix <- function (tb) {
    tb %>%
        mutate(row = factor(row, traits), col = factor(col, rev(traits))) %>%
        #mutate(correlation = ifelse(p_value < 0.05, correlation, NA)) %>%
        ggplot() +
        geom_tile(aes(x = row, y = col, fill = correlation)) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "maroon") +
        scale_x_discrete(position = "top", expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 0),
            axis.title = element_blank()
        ) +
        guides() +
        labs()

}

p <- treatments %>%
    calculate_cor() %>%
    plot_cor_matrix()

ggsave(paste0(folder_data, "temp/22-05-trait_correlation_matrix.png"), p, width = 8, height = 6)

# 6. dry weight by strain ----
p <- treatments %>%
    drop_na(strain) %>%
    ggplot() +
    geom_boxplot(aes(x = strain, y = dry_weight_mg, fill = strain_site_group), alpha = .6, outlier.size = -1, color = "black") +
    geom_point(aes(x = strain, y = dry_weight_mg, group = strain, color = strain), shape = 21, size = 2, stroke = 1, fill = NA,
               position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
    scale_color_manual(values = rep("black", 100)) +
    scale_fill_manual(values = rhizobia_site_colors, labels = c("High", "Low"), breaks = c("H", "L")) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.position = "none"
    ) +
    guides(color = "none") +
    labs(x = "", y = "above-ground dry weight (mg)")
ggsave(paste0(folder_data, "temp/22-06-weight_by_site.png"), p, width = 2, height = 4)


## Does rhizobia strain have effect on dry weight?
mod <- lmer(dry_weight_mg ~ strain + (1|plant) + (1|waterblock), data = treatments)
Anova(mod, type = 3) # rhizobia strain has effect on dry weight



# 7. Traits that respond to strain ----
tt_test <- treatments_long %>%
    nest(data = -trait) %>%
    mutate(fit = map(data, ~ lmer(value ~ strain + (1|plant) + (1|waterblock), data = .x)),
           ano = map(fit, ~Anova(.x, type = 3)),
           tidied = map(ano, tidy)) %>%
    unnest(tidied)

traits_sign <- tt_test %>%
    filter(term == "strain", p.value < 0.05) %>%
    pull(trait) # traits that are significant
tt <- treatments %>%
    select(id, strain, all_of(traits_sign)) %>%
    drop_na()
pcobj <- tt %>%
    select(-id, -strain) %>%
    prcomp(center = TRUE, scale. = TRUE)

p <- fviz_pca_ind(
    pcobj,
    label = "none",
    habillage = tt$strain,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
#    scale_color_manual(values = c(H = "#0C6291", L = "#BF4342"), labels = c("High", "Low"), breaks = c("H", "L"), name = "site") +
#    scale_shape_manual(values = c(H = 16, L = 17), labels = c("High", "Low"), breaks = c("H", "L"), name = "site") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(size = 5),
        legend.key.size = unit(3, "mm"),
        legend.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()
ggsave(paste0(folder_data, "temp/22-07-pca_strain.png"), p, width = 4, height = 4)


# 8. Table of all traits measured ----
features <- tibble(
    `Feature type` = c(rep("host yield", 2), rep("nodule", 2), rep("root architecture", 11)),
    `Feature` = traits
) %>%
    mutate(`Description` = c(
        "above-ground dry weight (mg)",
        "number of nodule",
        "root dry weight (mg)",
        "nodule dry weight (mg)",
        "number of root tips",
        "number of branching points",
        "total root length (px)",
        "branching frequency per px (1/px)",
        "area size (px^2)",
        "average diameter (px)",
        "median diameter (px)",
        "maximum diameter (px)",
        "estimated perimeter (px)",
        "estimated volume (px^3)",
        "estimated surface area (px^2)"
    )) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

ft <- features %>%
    flextable() %>%
    merge_at(i = 1:2, j = 2, part = "body") %>%
    merge_at(i = 3:4, j = 2, part = "body") %>%
    merge_at(i = 5:15, j = 2, part = "body") %>%
    width(j = 2, width = 2) %>%
    width(j = 3, width = 2) %>%
    width(j = 4, width = 4) %>%
    valign(j = 2, valign = "top", part = "all") %>%
    hline(i = 2, j = NULL, border = NULL, part = "body") %>%
    hline(i = 4, j = NULL, border = NULL, part = "body") %>%
    hline(i = 15, j = NULL, border = NULL, part = "body")

save_as_image(ft, paste0(folder_data, "temp/22-08-trait_table.png"), webshot = "webshot2")







# 9. correlation matrix ----

calculate_cor <- function (treatments) {
    temp <- treatments %>%
        select(all_of(traits)) %>%
        drop_na(all_of(traits)) %>%
        as.matrix() %>%
        rcorr(type = "spearman")

    tb_cor <- temp$r %>%
        as_tibble() %>%
        mutate(row = colnames(.)) %>%
        pivot_longer(-row, names_to = "col", values_to = "correlation")
    tb_p <- temp$P %>%
        as_tibble() %>%
        mutate(row = colnames(.)) %>%
        pivot_longer(-row, names_to = "col", values_to = "p_value")
    tb <- tb_cor %>%
        left_join(tb_p)
    return(tb)
}
plot_cor_matrix <- function (tb) {
    tb %>%
        mutate(row = factor(row, traits), col = factor(col, rev(traits))) %>%
        mutate(correlation = ifelse(p_value < 0.05, correlation, NA)) %>%
        ggplot() +
        geom_tile(aes(x = row, y = col, fill = correlation)) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "maroon") +
        scale_x_discrete(position = "top", expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 0),
            axis.title = element_blank()
        ) +
        guides() +
        labs()

}

p <- treatments %>%
    calculate_cor() %>%
    plot_cor_matrix()

ggsave(paste0(folder_data, "temp/22-09-trait_correlation_matrix.png"), p, width = 8, height = 6)


