#' This script computes the tradeoff between rhizobia growth vs. symbiosis traits

library(tidyverse)
library(janitor)
library(corrplot)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv"))

temp_clean <- function (tb) {
    tb %>%
        rename(gradient = population, population = site_group) %>%
        mutate(gradient = case_when(
            gradient == "VA" ~ "elevation",
            gradient == "PA" ~ "urbanization"
        ))
}
isolates <- temp_clean(isolates)
plants <- temp_clean(plants)

gst_r <- gts %>%
    group_by(exp_id) %>%
    mutate(temperature = factor(temperature, paste0(c(25, 30, 35, 40), "c"))) %>%
    filter(temperature == "30c") %>%
    arrange(temperature) %>%
    select(exp_id, temperature, r, lag, maxOD) %>%
    mutate(lag = -lag) %>% # change the sign of lag so greater lag = better growth
    pivot_wider(id_cols = exp_id, names_from = temperature, values_from = c(r, lag, maxOD)) %>%
    ungroup()

plants_r <- plants %>%
    group_by(gradient, population, exp_id, genome_id, exp_plant, exp_nitrogen) %>%
    select(shoot_biomass_mg, root_biomass_mg, shoot_height, primary_root_length, nodule_number, leaf_color, leaf_number) %>%
    summarize(across(everything(), function (x) mean(x, na.rm = T))) %>%
    filter(!exp_id == "control")

plants_rs <- list()
plants_rs[[1]] <- filter(plants_r, exp_plant == "lupulina", exp_nitrogen == "without nitrogen")
plants_rs[[2]] <- filter(plants_r, exp_plant == "sativa", exp_nitrogen == "without nitrogen")
plants_rs[[3]] <- filter(plants_r, exp_plant == "sativa", exp_nitrogen == "with nitrogen")

list_treatments <- c("lupulina w/o N", "sativa w/o N", "sativa w/ N")

for (i in 1:3) {
    isolates_trait <- isolates %>%
        filter(!genome_id %in% c("g2", "g3", "g15")) %>%
        #filter(gradient == "urbanization") %>%
        left_join(plants_rs[[i]]) %>%
        left_join(gst_r) %>%
        drop_na(exp_nitrogen) %>%
        ungroup()

    m <- isolates_trait[,-c(1:8)]
    m <- m[,(colSums(is.na(m)) == 0 & colSums(m) != 0)] # Remove traits with NA
    mm <- cor(m)
    p_mat <- cor.mtest(m, conf.level = 0.95)
    tt <- min(which(str_detect(names(m), "^r_"))) # the index of first growth trait; for plotting

    png(paste0(folder_data, "phenotypes/tradeoff/0", i, ".png"), width = 20, height = 20, units = "cm", res = 300)

    corrplot(
        mm, diag = T, title = paste0(list_treatments[i], ", n_strains = ", nrow(isolates_trait)), mar = c(0,0,2,0),
        p.mat = p_mat$p, insig = "label_sig", sig.level = c(0.001, 0.01, 0.05), pch.cex = 1,
        type = "full", tl.col = "black", tl.cex = 1, cl.ratio = 0.1, tl.srt = 45, tl.offset = .5
    ) %>%
        corrRect(index = c(1, tt, ncol(m)))


    dev.off()

}

