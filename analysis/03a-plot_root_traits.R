#' This script plots the root traits

library(tidyverse)
library(cowplot)
library(broom)
library(lme4) # for linear mixed-effect models
library(car) # companion to Applied Regression
library(DHARMa) # for checking assumptuons of lm4
library(broom.mixed) # for tidy up lme4
library(ggplotify) # for convert base plot into ggplot
library(factoextra) # for plotting pca eclipse
source(here::here("analysis/00-metadata.R"))

treatments <- read_csv(paste0(folder_data, "temp/03-treatments.csv"), show_col_types = F)
treatments_long <- read_csv(paste0(folder_data, "temp/03-treatments_long.csv"), show_col_types = F)
treatments_scaled <- read_csv(paste0(folder_data, "temp/03-treatments_scaled.csv"), show_col_types = F)
treatments_scaled_long <- read_csv(paste0(folder_data, "temp/03-treatments_scaled_long.csv"), show_col_types = F)

#
rhizobia_strains <- c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), unique(treatments$rhizobia))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
traits <- c("dry_weight", "nodule_number", "number_of_root_tips", "number_of_branch_points",
            "total_root_length_px", "branching_frequency_per_px", "network_area_px2",
            "average_diameter_px", "median_diameter_px", "maximum_diameter_px",
            "perimeter_px", "volume_px3", "surface_area_px2")

traits2 <- c(traits,
             paste0(rep(c("root_length_diameter_range_", "projected_area_diameter_range_", "surface_area_diameter_range_", "volume_diameter_range_"), each = 6),
                   rep(1:6, 4), rep(c("_px", "_px2", "_px2", "_px3"), each = 6)))

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
options(contrasts=c("contr.sum", "contr.poly"))


# 1. Traits by rhizobia location ----
# 1a. boxplot of traits by sites ----
p <- treatments_long %>%
    drop_na(value) %>%
    mutate(trait = factor(trait, traits)) %>%
    ggplot(aes(x = rhizobia_site, y = value, color = rhizobia_site)) +
    geom_boxplot(outlier.size = -1) +
    geom_jitter(width = .1, shape = 21) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y") +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/03a-01a-trait_site.png"), p, width = 10, height = 8)

##
treatments_scaled2 <- treatments_scaled %>%
    # Excluding the strains that do not nodulate
    filter(!rhizobia %in% c("H2M3R1", "L4M2R2"))
treatments_scaled_long2 <- treatments_scaled_long %>%
    # Excluding the strains that do not nodulate
    filter(!rhizobia %in% c("H2M3R1", "L4M2R2"))

# 1b. histogram of traits values ----
p <- treatments_long %>%
    mutate(trait = factor(trait, traits)) %>%
    drop_na(value) %>%
    group_by(trait) %>% mutate(n = n()) %>%
    ggplot() +
    geom_histogram(aes(x = value), color = "black", fill = "white") +
    geom_text(aes(x = Inf, y = Inf, label = paste0("n=", n)), hjust = 2, vjust = 2) +
    facet_wrap(~trait, scales = "free_x") +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/03a-01b-trait_value.png"), p, width = 15, height = 10)

# 1c. histogram of scaled traits values ----
p <- treatments_scaled_long2 %>%
    mutate(trait = factor(trait, traits)) %>%
    drop_na(value) %>%
    group_by(trait) %>% mutate(n = n()) %>%
    ggplot() +
    geom_histogram(aes(x = value), color = "black", fill = "white") +
    geom_text(aes(x = Inf, y = Inf, label = paste0("n=", n)), hjust = 2, vjust = 2) +
    facet_wrap(~trait, scales = "free_x") +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/03a-01c-trait_value_scaled.png"), p, width = 15, height = 10)

# 1d. check the LMM result ----
## Sanity check: does the biomass differ between plants with rhizobia from the two sites?
mod <- lmer(dry_weight ~ rhizobia_site + (1|rhizobia) + (1|waterblock) + (1|plant), data = treatments_scaled2)
Anova(mod, type = 3)
contrasts(as.factor(treatments_scaled2$rhizobia_site))

## LMM
traits_mod <- treatments_scaled_long2 %>%
    nest(data = -trait) %>%
    mutate(
        mod = map(data, ~ lmer(value ~ rhizobia_site + (1|rhizobia) + (1|waterblock) + (1|plant), data = .x)),
        ano = map(mod, ~ Anova(.x, type = 3)),
        tidied = map(ano, tidy),
        sr = map(mod, simulateResiduals)
        #tidied = map(mod, ~ tidy(.x, conf.int = T))
    ) %>%
    unnest(tidied) %>%
    filter(!str_detect(term, "Intercept")) %>%
    select(-data, -mod, -ano)
traits_mod

# 1e. check the model assumptions ----
p_list <- rep(list(NA), 13)
for (i in 1:13) p_list[[i]] <- as.ggplot(~plot(traits_mod$sr[[i]])) + theme(plot.background = element_rect(color = "black"))
p <- plot_grid(plotlist = p_list, nrow = 4, scale = .9, labels = traits, label_x = 0, hjust = 0) + theme(plot.background = element_rect(fill = "white"))
ggsave(paste0(folder_data, "temp/03a-01e-check_assumptions.png"), p, width = 30, height = 20)
# 1f. biomass as a function of rhizobia site in only M elevation plants ----
mod <- lmer(dry_weight ~ rhizobia_site + (1|rhizobia) + (1|waterblock) + (1|plant), data = treatments_scaled2 %>% filter(plant_site == "S"))
Anova(mod, type = 3)

tb <- treatments %>%
    drop_na(dry_weight) %>%
    filter(plant_site == "S")

tb_n <- tb %>%
    group_by(rhizobia_site) %>%
    count()

p <- tb %>%
    #group_by(rhizobia) %>% count
    #mutate(trait = factor(trait, traits)) %>%
    ggplot(aes(x = rhizobia_site, y = dry_weight, color = rhizobia_site)) +
    geom_boxplot(outlier.size = -1) +
    geom_jitter(width = .1, shape = 21) +
    geom_text(data = tb_n, aes(x = rhizobia_site, label = paste0("n=", n)), y = Inf, vjust = 1) +
    scale_color_manual(values = rhizobia_site_colors) +
    scale_y_continuous(expand = c(0.1, 0)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/03a-01f-trait_site.png"), p, width = 4, height = 3)




# 2. trait correlation ----
# 2a. dry_weight ~ nodule number, colored by trait site ----
p <- treatments %>%
    drop_na(rhizobia_site) %>%
    ggplot(aes(x = nodule_number, y = dry_weight, color = rhizobia_site)) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = rhizobia_site_colors) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/03a-02a-trait_site.png"), p, width = 4, height = 3)

## Does the rhizoibia site and nodule together have effect on biomass?
mod <- lmer(dry_weight ~ rhizobia_site + nodule_number + rhizobia_site:nodule_number +
                (1|rhizobia) + (1|waterblock) + (1|plant), data = treatments)
Anova(mod, type = 3) # nodule_number does, but not with rhizobia site

lmer(dry_weight ~ nodule_number + (1|rhizobia), data = treatments) %>%
    Anova(type = 3)


# 2b. correlation between all traits ----
library(GGally)
p <- treatments %>%
    drop_na(rhizobia_site) %>%
    #select(all_of(traits)) %>%
    ggpairs(
        columns = which(names(treatments) %in% traits),
        # columns = c(9,10,13),
        aes(color = rhizobia_site),
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", shape = 21)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5))
    )

ggsave(paste0(folder_data, "temp/03a-02b-trait_correlation.png"), p, width = 15, height = 15)

# 2c. correlation between all traits ----
p <- treatments_scaled %>%
    drop_na(rhizobia_site) %>%
    #select(all_of(traits)) %>%
    ggpairs(
        columns = which(names(treatments) %in% traits),
        # columns = c(9,10,13),
        aes(color = rhizobia_site),
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
        lower = list(continuous = wrap("smooth", shape = 21)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.5))
    )

ggsave(paste0(folder_data, "temp/03a-02c-trait_correlation_scaled.png"), p, width = 15, height = 15)

##
lmer(dry_weight ~ rhizobia_site + network_area_px2 + rhizobia_site:network_area_px2 +
         (1|rhizobia) + (1|waterblock) + (1|plant), data = treatments) %>%
    Anova(type = 3)

# 2d. Correlation matrix with significant level ----
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
    scale_x_discrete(position = "top",) +
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

ggsave(paste0(folder_data, "temp/03a-02d-trait_correlation_matrix.png"), p, width = 10, height = 8)
# 2e. correlation matrix by strain ----
p_list <- rep(list(NA), 6)
for (i in 1:6) {
    p_list[[i]] <- treatments %>%
        filter(rhizobia == rhizobia_strains[i]) %>%
        calculate_cor() %>%
        plot_cor_matrix()
}

a <- table(treatments$rhizobia)
p <- plot_grid(plotlist = p_list, nrow = 2, labels = paste0(names(a), " (n=", a, ")"), hjust = 0)
ggsave(paste0(folder_data, "temp/03a-02e-trait_correlation_matrix_strain.png"), p, width = 30, height = 16)


# 3. Traits by strain ----
# 3a. histogram of traits by strains
p <- treatments_scaled_long %>%
    drop_na(value) %>%
    mutate(trait = factor(trait, traits)) %>%
    ggplot(aes(x = rhizobia, y = value, color = rhizobia_site)) +
    geom_boxplot() +
    geom_jitter(width = .1, shape = 21) +
    scale_color_manual(values = rhizobia_site_colors) +
    facet_wrap(~trait, scales = "free_y") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/03a-03a-trait_strain.png"), p, width = 10, height = 8)

# 3b. LMM result ----
## Sanity check: does the biomass differ between plants with rhizobia from the two sites?
mod <- lmer(dry_weight ~ rhizobia + (1|waterblock) + (1|plant), data = treatments_scaled)
Anova(mod, type = 3)

## LMM
traits_mod <- treatments_scaled_long %>%
    nest(data = -trait) %>%
    mutate(
        mod = map(data, ~ lmer(value ~ rhizobia + (1|waterblock) + (1|plant), data = .x)),
        ano = map(mod, ~ Anova(.x, type = 3)),
        tidied = map(ano, tidy),
        sr = map(mod, simulateResiduals)
    ) %>%
    unnest(tidied) %>%
    filter(!str_detect(term, "Intercept")) %>%
    select(-data, -mod, -ano)
traits_mod
traits_mod %>%
    filter(p.value < 0.05)

p_list <- rep(list(NA), 13)
for (i in 1:13) p_list[[i]] <- as.ggplot(~plot(traits_mod$sr[[i]])) + theme(plot.background = element_rect(color = "black"))
p <- plot_grid(plotlist = p_list, nrow = 4, scale = .9, labels = traits, label_x = 0, hjust = 0) + theme(plot.background = element_rect(fill = "white"))
ggsave(paste0(folder_data, "temp/03a-03b-check_assumptions.png"), p, width = 30, height = 20)

#
options(contrasts = c("contr.treatment", "contr.treatment"))

# 3c. LM of traits ----
## Sanity check
lm(dry_weight ~ rhizobia, data = treatments_scaled) %>%
    summary()

lm(nodule_number ~ rhizobia, data = treatments_scaled) %>%
    summary()

## OLM
traits_mod <- treatments_scaled_long %>%
    mutate(rhizobia = factor(rhizobia, c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2"))) %>%
    nest(data = -trait) %>%
    mutate(
        mod = map(data, ~ lm(value ~ rhizobia, data = .x)),
        tidied = map(mod, tidy)
    ) %>%
    unnest(tidied) %>%
    select(-data, -mod)

## OLM significance
traits_mod %>%
    mutate(significance = case_when(
        p.value > 0.01 & p.value < 0.05 ~ "*",
        p.value > 0.001 & p.value < 0.01 ~ "**",
        p.value < 0.001 ~ "***",
        TRUE ~ "n.s."
    )) %>%
    select(trait, term, significance) %>%
    group_by(trait) %>%
    pivot_wider(names_from = term, values_from = significance)


# 4. PCA ----
# 4a. pca ----
compute_pca_coord <- function (pcobj) {
    #' The function below comes from the source code of ggbiplot
    #' to replace the use of ggbiplot
    choices = 1:2 # which PCs to plot
    scale = 1
    obs.scale = 1 - scale
    var.scale = scale

    # Recover the SVD
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
    # Scores
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
    # Directions
    v <- sweep(v, 2, d^var.scale, FUN='*')
    df.v <- as.data.frame(v[, choices])

    names(df.u) <- c('xvar', 'yvar')
    names(df.v) <- names(df.u)
    df.u <- df.u * nobs.factor

    # Variable Names
    df.v$varname <- rownames(v)
    # Variables for text label placement
    df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
    df.v$hjust = with(df.v, (1 - 1.5 * sign(xvar)) / 2)

    # Change the labels for the axes
    # Append the proportion of explained variance to the axis labels
    u.axis.labs <- paste('standardized PC', choices, sep='')
    u.axis.labs <- paste(u.axis.labs, sprintf('(%0.1f%% explained var.)', 100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))

    return(list(df.v = tibble(df.v), # Variables
                df.u = tibble(df.u), # Score
                u.axis.labs = u.axis.labs))
}
pcobj <- treatments %>%
    select(id, all_of(traits)) %>%
    #select(id, all_of(c("dry_weight", "nodule_number"))) %>%
    drop_na() %>%
    select(-id) %>%
    prcomp(center = TRUE, scale. = TRUE)
pca_coord <- compute_pca_coord(pcobj)

p <- pca_coord$df.u %>%
    bind_cols(drop_na(treatments, all_of(traits))) %>%
    #arrange(ColorLabel) %>%
    ggplot() +
    # Draw directions of axis
    geom_segment(data = pca_coord$df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 'picas')), color = scales::muted('red')) +
    # Draw scores
    geom_point(aes(x = xvar, y = yvar, color = rhizobia_site), shape = 21,
                   #color = ColorLabel, fill = FillLabel, shape = ShapeLabel, alpha = AlphaLabel),
               size = 2, stroke = .8) +
    # Label the variable axes
    geom_text(data = pca_coord$df.v, aes(label = varname, x = xvar, y = yvar, angle = angle, hjust = hjust), color = 'blue', size = 3) +
    # scale_color_manual(values = color_names, name = "", label = names(color_names)) +
    # scale_fill_manual(values = fill_names, name = "", label = names(color_names)) +
    # scale_shape_manual(values = shape_names, name = "", label = names(color_names)) +
    # scale_alpha_manual(values = alpha_names, name = "", label = names(color_names)) +
    theme_classic() +
    labs(x = pca_coord$u.axis.labs[1], y = pca_coord$u.axis.labs[2])

ggsave(paste0(folder_data, "temp/03a-04a-pca.png"), p, width = 8, height = 7)

# 4b. four strains that nodulate
nodulating_strains <- c("H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1")
pcobj <- treatments %>%
    filter(rhizobia %in% nodulating_strains) %>%
    select(id, all_of(traits)) %>%
    #select(id, all_of(c("dry_weight", "nodule_number"))) %>%
    drop_na() %>%
    select(-id) %>%
    prcomp(center = TRUE, scale. = TRUE)
pca_coord <- compute_pca_coord(pcobj)

p <- pca_coord$df.u %>%
    bind_cols(drop_na(filter(treatments, rhizobia %in% nodulating_strains), all_of(traits))) %>%
    ggplot() +
    # Draw directions of axis
    geom_segment(data = pca_coord$df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 'picas')), color = scales::muted('red')) +
    # Draw scores
    geom_point(aes(x = xvar, y = yvar, color = rhizobia), shape = 21,
               size = 2, stroke = .8) +
    # Label the variable axes
    geom_text(data = pca_coord$df.v, aes(label = varname, x = xvar, y = yvar, angle = angle, hjust = hjust), color = 'blue', size = 3) +
    theme_classic() +
    labs(x = pca_coord$u.axis.labs[1], y = pca_coord$u.axis.labs[2])

ggsave(paste0(folder_data, "temp/03a-04b-pca_nodulating_strains.png"), p, width = 8, height = 7)
p <- fviz_pca_ind(pcobj, habillage = drop_na(filter(treatments, rhizobia %in% nodulating_strains), all_of(traits))$rhizobia, addEllipses = T)
ggsave(paste0(folder_data, "temp/03a-04b-pca_nodulating_strains_ec.png"), p, width = 8, height = 7)

# 4b. six strains omitting nodule number
pcobj <- treatments %>%
    select(id, all_of(traits[-2])) %>%
    drop_na() %>%
    select(-id) %>%
    prcomp(center = TRUE, scale. = TRUE)
pca_coord <- compute_pca_coord(pcobj)

p <- pca_coord$df.u %>%
    bind_cols(drop_na(treatments, all_of(traits[-2]))) %>%
    ggplot() +
    # Draw directions of axis
    geom_segment(data = pca_coord$df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 'picas')), color = scales::muted('red')) +
    # Draw scores
    geom_point(aes(x = xvar, y = yvar, color = rhizobia), shape = 21,
               size = 2, stroke = .8) +
    # Label the variable axes
    geom_text(data = pca_coord$df.v, aes(label = varname, x = xvar, y = yvar, angle = angle, hjust = hjust), color = 'blue', size = 3) +
    theme_classic() +
    labs(x = pca_coord$u.axis.labs[1], y = pca_coord$u.axis.labs[2])

ggsave(paste0(folder_data, "temp/03a-04c-pca_strains.png"), p, width = 8, height = 7)

p <- fviz_pca_ind(pcobj, habillage = drop_na(treatments, all_of(traits[-2]))$rhizobia, addEllipses = T)
ggsave(paste0(folder_data, "temp/03a-04c-pca_strains_ec.png"), p, width = 8, height = 7)
#


#
install.packages("factoextra")
library(factoextra)
iris.pca <- prcomp(iris[, -5],
                   scale = TRUE)

summary(iris.pca)

fviz_pca_ind(iris.pca,
             habillage=iris$Species,
             addEllipses=TRUE)















