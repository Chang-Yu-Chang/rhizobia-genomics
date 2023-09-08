#' This script plots the root trait

library(tidyverse)
library(cowplot)
library(broom)
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

# plant traits
treatments <- read_csv(paste0(folder_data, "temp/11-treatments.csv"), show_col_types = F)

# Clean up data. Use only medium elevation plants
treatments_M <- treatments %>%
    mutate(strain_site_group = ifelse(is.na(strain_site_group), "control", strain_site_group),
           strain_site = ifelse(is.na(strain_site), "control", strain_site),
           strain = ifelse(is.na(strain), "control", strain)) %>%
    filter(plant_site_group == "S") %>%
    mutate(strain_site_group = factor(strain_site_group, c("H", "L", "control")))

# 1. Is high growth rate associated with high root biomass at the strain level?
#rhizobia_strains <- c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")
gc1 <- gc.prm %>% filter(strain %in% rhizobia_strains)
tm1 <- treatments_M %>% filter(strain != "control")


n_bs <- 1000 # number of bootstraps
n_spb <- 10 # number of resamples per strain per bootstrap
list_cor <- tibble(r0 = rep(NA, n_bs), r1 = rep(NA, n_bs))
n_bsp <- n_spb * length(rhizobia_strains) # number of sampled pairs in one bootstrap
list_pair <- rep(list(tibble(gi0 = rep(NA, n_bsp), ei0 = rep(NA, n_bsp), gi1 = rep(NA, n_bsp), ei1 = rep(NA, n_bsp))), n_bs)

for (i in 1:n_bs) {
    set.seed(i)
    # Resample across strains. This is the randomized control
    gi0 <- sample(gc1$r[!is.na(gc1$r)], n_spb * length(rhizobia_strains), replace = T)
    ei0 <- sample(tm1$dry_weight_mg[!is.na(tm1$dry_weight_mg)], n_spb * length(rhizobia_strains), replace = T)
    list_cor$r0[i] <- cor(gi0, ei0)


    # Resample within strains
    gi1 <- rep(list(NA), length(rhizobia_strains))
    ei1 <- rep(list(NA), length(rhizobia_strains))
    for (j in 1:length(rhizobia_strains)) {
        #' Edit here to make it a function so that we can use different traits as arguments
        gj <- gc1$r[gc1$strain == rhizobia_strains[j]] # growth trait.
        ej <- tm1$dry_weight_mg[tm1$strain == rhizobia_strains[j]] # extended phenotype
        gi1[[j]] <- sample(gj[!is.na(gj)], n_spb, replace = T)
        ei1[[j]] <- sample(ej[!is.na(ej)], n_spb, replace = T)
    }

    gi1 <- unlist(gi1)
    ei1 <- unlist(ei1)
    list_cor$r1[i] <- cor(gi1, ei1)

    # Save the bootstrapped samples
    list_pair[[i]]$gi0 <- gi0
    list_pair[[i]]$ei0 <- ei0
    list_pair[[i]]$gi1 <- gi1
    list_pair[[i]]$ei1 <- ei1

    cat(" ", i)

}

list_cor <- list_cor %>%
    mutate(bootstrap = 1:n()) %>%
    mutate(difference = ifelse(r1 > r0, "r1 is greater", "r1 is less"))

boot_pair <- bind_rows(list_pair, .id = "bootstrap")
write_csv(boot_pair, paste0(folder_data, "temp/12-01-raw_growth_rate_vs_biomass.csv"))


table(list_cor$difference)
mean(list_cor$r0) # 0.009594945
mean(list_cor$r1) # 0.4390422
# r0 is randomized_across_strains
# r1 is randomized_within_strains
p <- list_cor %>%
    pivot_longer(cols = c(-bootstrap, -difference), names_to = "treatment") %>%
    #mutate(treatment = ifelse(treatment == "r0", "randomized_across_strains", "randomized_within_strains")) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = treatment, y = value), shape = 21) +
    geom_segment(data = list_cor, x = "r0", xend = "r1", aes(y = r0, yend = r1, color = difference, group = bootstrap), linewidth = 0.1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 3, "Set1")) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "")

ggsave(paste0(folder_data, "temp/12-01-growth_rate_vs_biomass.png"), p, width = 4, height = 3)


# 2. growth rate vs. nodule number
list_cor <- tibble(r0 = rep(NA, n_bs), r1 = rep(NA, n_bs))
for (i in 1:n_bs) {
    set.seed(i)
    # Resample across strains. This is the randomized control
    gi0 <- sample(gc1$r[!is.na(gc1$r)], n_spb * length(rhizobia_strains), replace = T)
    ei0 <- sample(tm1$nodule_number[!is.na(tm1$nodule_number)], n_spb * length(rhizobia_strains), replace = T)
    list_cor$r0[i] <- cor(gi0, ei0)


    # Resample within strains
    gi1 <- rep(list(NA), length(rhizobia_strains))
    ei1 <- rep(list(NA), length(rhizobia_strains))
    for (j in 1:length(rhizobia_strains)) {
        gj <- gc1$r[gc1$strain == rhizobia_strains[j]] # growth trait.
        ej <- tm1$nodule_number[tm1$strain == rhizobia_strains[j]] # extended phenotype
        gi1[[j]] <- sample(gj[!is.na(gj)], n_spb, replace = T)
        ei1[[j]] <- sample(ej[!is.na(ej)], n_spb, replace = T)
    }

    gi1 <- unlist(gi1)
    ei1 <- unlist(ei1)
    list_cor$r1[i] <- cor(gi1, ei1)


    cat(" ", i)

}

list_cor <- list_cor %>%
    mutate(bootstrap = 1:n()) %>%
    mutate(difference = ifelse(r1 > r0, "r1 is greater", "r1 is less"))

# r0 is randomized_across_strains
# r1 is randomized_within_strains
p <- list_cor %>%
    pivot_longer(cols = c(-bootstrap, -difference), names_to = "treatment") %>%
    #mutate(treatment = ifelse(treatment == "r0", "randomized_across_strains", "randomized_within_strains")) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = treatment, y = value), shape = 21) +
    geom_segment(data = list_cor, x = "r0", xend = "r1", aes(y = r0, yend = r1, color = difference, group = bootstrap), linewidth = 0.1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 3, "Set1"), breaks = c("r1 is greater", "r1 is less")) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "")

ggsave(paste0(folder_data, "temp/12-02-growth_rate_vs_nodule.png"), p, width = 4, height = 3)


















