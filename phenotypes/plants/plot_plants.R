#' This script plots the trait data of plant experiments

library(tidyverse)
library(cowplot)
library(ggsci)
library(vegan)
# library(lme4) # for linear mixed-effect models
# library(car) # companion to Applied Regression
# library(broom.mixed) # for tidy up lme4
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_phenotypes, "plants/plants.csv"))

# Summary
nrow(plants) # 822 plants
table(plants$exp_plant) # 251 lupulina, 571 sativa
length(unique(plants$exp_id)) # 26 rhizobia strains used in the plant experiments

# 1. plot all plant traits ----
set.seed(1)
p <- plants %>%
    filter(population != "control") %>%
    filter(exp_plant == "sativa", exp_nitrogen == "without nitrogen") %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    #filter(str_detect(trait, "shoot|leaf|nodule")) %>%
    ggplot(aes(x = population, y = value, fill = population)) +
    geom_boxplot(outlier.shape = -1, alpha = 0.5) +
    geom_jitter(size = .5, shape = 21, width = .2) +
    scale_fill_manual(values = population_colors) +
    coord_cartesian(clip = "off") +
    facet_wrap(gradient~trait, scales = "free", ncol = 7) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)
       # strip.background = element_rect(color = NA, fill = "grey90")
    ) +
    guides(fill = "none") +
    labs(y = "")

ggsave(paste0(folder_phenotypes, "plants/01-sativa_traits.png"), p, width = 14, height = 8)

p <- plants %>%
    filter(population != "control") %>%
    filter(exp_plant == "sativa", gradient == "elevation") %>%
    mutate(exp_nitrogen = factor(exp_nitrogen, c("without nitrogen", "with nitrogen"))) %>%
    select(-nodule_shape, -nodule_size, -nodule_color, -exp_labgroup) %>%
    group_by(gradient, population, exp_plant) %>%
    filter(nodule_number <100) %>%
    pivot_longer(cols = -c(1:11), names_to = "trait", values_drop_na = T) %>%
    #filter(str_detect(trait, "shoot|leaf|nodule")) %>%
    ggplot(aes(x = population, y = value, fill = population)) +
    geom_boxplot(outlier.shape = -1, alpha = 0.5) +
    geom_jitter(size = .5, shape = 21, width = .2) +
    scale_fill_manual(values = population_colors) +
    coord_cartesian(clip = "off") +
    facet_wrap(exp_nitrogen~trait, scales = "free", ncol = 7) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)
        # strip.background = element_rect(color = NA, fill = "grey90")
    ) +
    guides(fill = "none") +
    labs(y = "")

ggsave(paste0(folder_phenotypes, "plants/01-sativa_nitrogen.png"), p, width = 14, height = 8)

# 2. Nodule vs. plant biomass ----
p <- plants %>%
    filter(exp_plant == "sativa") %>%
    filter(population != "control") %>%
    filter(nodule_number <100) %>%
    #select(nodule_number, shoot_biomass_mg)
    #filter(str_detect(trait, "shoot|leaf|nodule")) %>%
    ggplot() +
    geom_point(aes(x = nodule_number, y = shoot_height, color = exp_nitrogen, shape = exp_nitrogen), size = 1, stroke = 1) +
    scale_shape_manual(values = c(`with nitrogen` = 19, `without nitrogen` = 21), labels = c("N+", "N-")) +
    scale_color_aaas(labels = c("N+", "N-")) +
    #geom_point(aes(x = nodule_number, y = shoot_biomass_mg), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
    #scale_x_continuous(breaks = seq(0, 50, 10)) +
    #coord_equal() +
    #facet_grid(trait~population, scales = "free", space = "free_x") +
    theme_bw() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_phenotypes, "plants/02-nodule_vs_height.png"), p, width = 6, height = 5)

# 3. PCAs ----
tbs <- tibble(
    gradient = c("elevation", "urbanization", "elevation", "urbanization", "elevation"),
    exp_plant = c("lupulina", "lupulina", "sativa", "sativa", "sativa"),
    exp_nitrogen = c("N-", "N-", "N-", "N-", "N+")
) %>%
    mutate(treatment = paste(gradient, exp_plant, exp_nitrogen, sep = "\t")) %>%
    mutate(plants_data = list(
        # Lupulinas. First plant exp
        plants1 = plants %>%
            filter(exp_plant == "lupulina", exp_id != "control", gradient == "elevation") %>%
            select(population, exp_id, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
            drop_na(),
        plants2 =  plants %>%
            filter(exp_plant == "lupulina", exp_id != "control", gradient == "urbanization") %>%
            select(population, exp_id, shoot_biomass_mg, root_biomass_mg, nodule_number) %>%
            drop_na(),
        # sativas. Second plant exp
        plants3 = plants %>%
            filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", gradient == "elevation") %>%
            select(population, exp_id, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>%
            drop_na(),
        plants4 = plants %>%
            filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "without nitrogen", gradient == "urbanization") %>%
            select(population, exp_id, shoot_height, nodule_number, leaf_number, leaf_color, lateral_root_number, longest_lateral_root_length) %>%
            drop_na(),
        plants5 = plants %>%
            filter(exp_plant == "sativa", exp_id != "control", exp_nitrogen == "with nitrogen", gradient == "elevation") %>%
            select(population, exp_id, shoot_height, nodule_number, longest_petiole_length, leaf_number, leaf_color) %>%
            drop_na()
    ))


do_pca <- function(x) {
    prcomp(x[,-c(1,2)], scale. = TRUE)
}
get_pcs <- function (plants_subset, pca_result) {
    #' Extract the PCs from the pca object
    as_tibble(`[[`(pca_result, "x")) %>%
        mutate(population = plants_subset$population, exp_id = plants_subset$exp_id) %>%
        left_join(distinct(isolates, gradient, population))
}
get_pcvar <- function (pca_result) {
    #' Get the variance explained by the top two PCs
    summary(pca_result)$importance[2, ] %>% round(3) * 100
}
plot_pca <- function (pcsi, pca_result) {
    dm <- vegdist(select(pcsi, starts_with("PC")), method = "euclidean")
    # strata by exp_id
    #mod <- with(pcsi, adonis2(dm ~ population, data = pcsi, permutations = 10000, strata = exp_id))
    mod <- adonis2(dm ~ population, data = pcsi, permutations = 10000)
    pcsi %>%
        ggplot() +
        geom_point(aes(x = PC1, y = PC2, color = population), shape = 21, stroke = 1, size = 2) +
        stat_ellipse(aes(x = PC1, y = PC2, fill = population), geom = "polygon", type = "norm", level = 0.95, alpha = .2) +
        annotate("text", x = Inf, y = Inf, label = paste0("N=", nrow(pcsi)), hjust = 1.1, vjust = 1.1) +
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, label = paste0("p=", round(mod[1,5], 5))) +
        geom_vline(xintercept = 0, color = "grey10", linetype = 2) +
        geom_hline(yintercept = 0, color = "grey10", linetype = 2) +
        scale_color_manual(values = population_colors) +
        scale_fill_manual(values = population_colors, name = "population") +
        scale_x_continuous(breaks = seq(-8,8,2)) +
        scale_y_continuous(breaks = seq(-8,8,2)) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey10", fill = NA),
            legend.background = element_blank()
            #plot.background = element_blank(),
            #plot.margin = margin(0,0,0,0, "mm")
        ) +
        guides(color = "none") +
        labs(x = paste0("PC1 (", get_pcvar(pca_result)[1], "%)"), y = paste0("PC2 (", get_pcvar(pca_result)[2], "%)"))

}

tbs_pca <- tbs %>%
    mutate(
        pca_results = map(plants_data, do_pca),
        pcs = map2(plants_data, pca_results, get_pcs),
        p_pca = map2(pcs, pca_results, plot_pca)
    )
p_pcas <- tbs_pca$p_pca

p <- plot_grid(
    p_pcas[[3]] + ggtitle("elevation sativa N-") + theme(legend.position = "none"),
    p_pcas[[5]] + ggtitle("elevation sativa N+") + theme(legend.position = "bottom", legend.title = element_blank()),
    nrow = 1, scale = .95, align = "vh", axis = "lr") +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_phenotypes, "plants/03-pca_nitrogen.png"), p, width = 9, height = 5)


p <- plot_grid(
    p_pcas[[1]] + ggtitle("elevation lupulina N-") + theme(legend.position = "none"),
    p_pcas[[2]] + ggtitle("urbanization lupulina N-") + theme(legend.position = "none"),
    p_pcas[[3]] + ggtitle("elevation sativa N-") + theme(legend.position = "bottom", legend.title = element_blank()),
    p_pcas[[4]] + ggtitle("urbanization sativa N-") + theme(legend.position = "bottom", legend.title = element_blank()),
    nrow = 2, scale = .95, align = "vh", axis = "lr") +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_phenotypes, "plants/03-pcas.png"), p, width = 9, height = 9)


if (F) {

plants %>%
    arrange(site_group) %>%
    distinct(id, site_group, exp_id) %>%
    group_by(site_group, exp_id) %>%
    count()

p <- plants_long %>%
    ggplot() +
    geom_boxplot(aes(x = exp_id, y = value, color = site_group), outlier.shape = NA) +
    geom_jitter(aes(x = exp_id, y = value, color = site_group), width = 0.2, height = 0, shape = 21) +
    facet_grid(trait~site_group, scales = "free", space = "free_x") +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    guides(color = "none") +
    labs()

ggsave(paste0(folder_phenotypes, "symbiosis/01-plant_traits.png"), p, width = 10, height = 10)

}



# deprecated ----
if (F) {

# 2. Plot the nodule vs plant biomass
p <- plants_inoculated %>%
    ggplot() +
    geom_smooth(aes(x = nodule_count, y = shoot_biomass_mg), method = "lm")+
    geom_point(aes(x = nodule_count, y = shoot_biomass_mg), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    coord_equal() +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey95")
    ) +
    guides() +
    labs(x = "# of nodules", y = "shoot biomass (mg)")

ggsave(paste0(folder_phenotypes, "symbiosis/02-nodule_vs_shoot.png"), p, width = 5, height = 4)

# Stats
lm(shoot_biomass_mg ~ nodule_count, data = plants) %>% broom::tidy()

# 3. Plot the pairwise comparison of the three symbiosis traits
plot_dots <- function(tb_wide, trait1, trait2) {
    tb_wide %>%
        ggplot() +
        geom_smooth(aes(x = {{trait1}}, y = {{trait2}}), method = "lm")+
        geom_point(aes(x = {{trait1}}, y = {{trait2}}), shape = 21, size = 2, alpha = 0.5, stroke = 1)+
        theme_classic() +
        theme(
            panel.grid.major = element_line(color = "grey95")
        ) +
        guides()
}
p1 <- plants %>%
    filter(site_group != "control") %>%
    plot_dots(nodule_count, shoot_biomass_mg) +
    labs(x = "# of nodules", y = "shoot biomass (mg)")
p2 <- plants %>%
    filter(site_group != "control") %>%
    plot_dots(nodule_count, root_biomass_mg) +
    labs(x = "# of nodules", y = "root biomass (mg)")
p3 <- plants %>%
    filter(site_group != "control") %>%
    plot_dots(root_biomass_mg, shoot_biomass_mg) +
    labs(x = "root biomass (mg)", y = "shoot biomass (mg)")
p <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb", scale = 0.95, labels = LETTERS[1:3]) +
    theme(plot.background = element_rect(color = NA, fill = "white"))

ggsave(paste0(folder_phenotypes, "symbiosis/03-symbiosis_trait_pairs.png"), p, width = 9, height = 3)

# Are the symbiosis traits mutually correlated?
cor.test(plants$nodule_count, plants$shoot_biomass_mg) %>% broom::tidy()
cor.test(plants$nodule_count, plants$root_biomass_mg) %>% broom::tidy()
cor.test(plants$root_biomass_mg, plants$shoot_biomass_mg) %>% broom::tidy()

# 4. Plot the symbiosis traits comparing the strains
#' Plot raw data and sem
compute_trait_mean <- function (plants_long, tra = "shoot_biomass_mg", pop = "VA") {
    pl <- plants_long %>%
        filter(trait == tra) %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, exp_id, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = exp_id, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = exp_id, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = exp_id, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "shoot biomass (mg)")
}

t1 <- compute_trait_mean(plants_long, tra = "shoot_biomass_mg", pop = "VA")
t2 <- compute_trait_mean(plants_long, tra = "shoot_biomass_mg", pop = "PA")
p1 <- plot_dots(t1$pl, t1$plm)
p2 <- plot_dots(t2$pl, t2$plm)
p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("", ""))
ggsave(paste0(folder_phenotypes, "symbiosis/04-trait_strain.png"), p, width = 5, height = 3)

# # Stats
# plants_test <- filter(plants, population == "VA")
# ## Does rhizobia strain have effect on shoot biomass in VA?
# mod <- lmer(shoot_biomass_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yup
# ## Does rhizobia strain have effect on nodule number in VA?
# mod <- lmer(nodule_count ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yup
# ## Does rhizobia strain have effect on root biomass in VA?
# mod <- lmer(root_biomass_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # nope
#
# plants_test <- filter(plants, population == "PA")
# ## Does rhizobia strain have effect on shoot biomass in PA?
# mod <- lmer(shoot_biomass_mg ~ exp_id + (1|site_group) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes
# ## Does rhizobia strain have effect on nodule number in VA?
# mod <- lmer(nodule_count ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes
# ## Does rhizobia strain have effect on root biomass in VA?
# mod <- lmer(root_biomass_mg ~ exp_id + (1|site_group) + (1|plant) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # nope


# 5. Plot the symbiosis traits comparing the two populations
compute_trait_mean2 <- function (plants_long, tra = "shoot_biomass_mg", pop = "VA") {
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
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
        scale_fill_manual(values = site_group_colors) +
        facet_wrap(.~site_group, scales = "free_x", nrow = 1) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "shoot biomass (mg)")
}

t1 <- compute_trait_mean2(plants_long,tra = "shoot_biomass_mg", pop = "VA")
t2 <- compute_trait_mean2(plants_long,tra = "shoot_biomass_mg", pop = "PA")
p1 <- plot_dots2(t1$pl, t1$plm)
p2 <- plot_dots2(t2$pl, t2$plm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_phenotypes, "symbiosis/05-trait_population.png"), p, width = 5, height = 3)

# # Does rhizobia sites have effect on shoot biomass?
# plants_test <- filter(plants, population == "VA")
# mod <- lmer(shoot_biomass_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # no
# mod <- lmer(nodule_count ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # no
# mod <- lmer(root_biomass_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # no
#
# # Does rhizobia sites have effect on shoot biomass?
# plants_test <- filter(plants, population == "PA")
# mod <- lmer(shoot_biomass_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes
# mod <- lmer(nodule_count ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes
# mod <- lmer(root_biomass_mg ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # yes

# 6. Plot the PCA for the four sites
pw1 <- plants %>%
    drop_na(shoot_biomass_mg, nodule_count, root_biomass_mg) %>%
    filter(population == "VA", exp_id != "control")
pcobj1 <- pw1 %>%
    select(shoot_biomass_mg, nodule_count, root_biomass_mg) %>%
    prcomp(center = TRUE, scale. = TRUE)

p1 <- fviz_pca_ind(
    pcobj1,
    label = "none",
    habillage = pw1$site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "top",
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()

pw2 <- plants %>%
    drop_na(shoot_biomass_mg, nodule_count, root_biomass_mg) %>%
    filter(population == "PA", exp_id != "control")
pcobj2 <- pw2 %>%
    select(shoot_biomass_mg, nodule_count, root_biomass_mg) %>%
    prcomp(center = TRUE, scale. = TRUE)

p2 <- fviz_pca_ind(
    pcobj2,
    label = "none",
    habillage = pw2$site_group,
    addEllipses = TRUE, ellipse.level = 0.95, ellipse.alpha = 0
) +
    scale_color_manual(values = site_group_colors) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "top",
        legend.background = element_rect(fill = NA, color = NA),
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_blank()
    ) +
    guides(fill = "none") +
    labs()

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"), scale = 0.9) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(paste0(folder_phenotypes, "symbiosis/06-trait_pca.png"), p, width = 9, height = 5)

# # Does rhizobia sites have effect on PC1 of symbiosis in VA?
# plants_test <- mutate(pw1, pc1 = pcobj1$x[,1]) %>% left_join(plants)
# mod <- lmer(pc1 ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3) # No difference between high and low elevation
#
# # Does rhizobia sites have effect on PC1 of symbiosis in VA?
# plants_test <- mutate(pw2, pc1 = pcobj2$x[,1]) %>% left_join(plants)
# mod <- lmer(pc1 ~ site_group + (1|exp_id) + (1|waterblock), data = plants_test)
# Anova(mod, type = 3)  # No differnece between urban and suburban
# #pcobj2$sdev^2 / sum(pcobj2$sdev^2)


# 7. Plot all three traits in dot plots by strain
compute_trait_mean3 <- function (plants_long, pop = "VA") {
    pl <- plants_long %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, exp_id, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots3 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = exp_id, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = exp_id, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = exp_id, ymin = mean-sem, ymax = mean+sem), width = 0.5) +
        scale_fill_manual(values = site_group_colors) +
        facet_grid(trait~site_group, scales = "free") +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "")
}

t1 <- compute_trait_mean3(plants_long, pop = "VA")
t2 <- compute_trait_mean3(plants_long, pop = "PA")
p1 <- plot_dots3(t1$pl, t1$plm)
p2 <- plot_dots3(t2$pl, t2$plm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_phenotypes, "symbiosis/07-trait_strain.png"), p, width = 8, height = 9)


# 6. Plot all three traits in violin plot
compute_trait_mean4 <- function (plants_long, pop = "VA") {
    pl <- plants_long %>%
        filter(exp_id != "control") %>%
        filter(population == pop)

    plm <- pl %>%
        group_by(population, site_group, trait) %>%
        summarize(mean = mean(value), sem = sd(value)/sqrt(n()))
    return(list(pl = pl, plm = plm))
}
plot_dots4 <- function (pl, plm) {
    set.seed(1)
    pl %>%
        ggplot() +
        geom_rect(aes(fill = site_group), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
        geom_jitter(aes(x = site_group, y = value), alpha = 0.3, shape = 16, color = "black", height = 0, width = 0.1) +
        geom_point(data = plm, aes(x = site_group, y = mean), size = 2, shape = 21, fill = NA, color = "black") +
        geom_errorbar(data = plm, aes(x = site_group, ymin = mean-sem, ymax = mean+sem), width = 0.2) +
        scale_fill_manual(values = site_group_colors) +
        facet_grid(trait~site_group, scales = "free") +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        guides(fill = "none") +
        labs(x = " ", y = "")
}

t1 <- compute_trait_mean4(plants_long, "VA")
t2 <- compute_trait_mean4(plants_long, "PA")
p1 <- plot_dots4(t1$pl, t1$plm)
p2 <- plot_dots4(t2$pl, t2$plm)

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", labels = c("VA", "PA"))
ggsave(paste0(folder_phenotypes, "symbiosis/08-trait_population.png"), p, width = 8, height = 9)

}
