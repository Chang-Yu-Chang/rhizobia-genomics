#' This script computes the tradeoff between rhizobia growth vs. symbiosis traits

library(tidyverse)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
plants <- read_csv(paste0(folder_data, "phenotypes/plants/plants.csv"))
gts <- read_csv(paste0(folder_data, "phenotypes/growth/gts.csv"))

# Traits
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
    #select(shoot_biomass_mg, root_biomass_mg, shoot_height, nodule_number, leaf_color, leaf_number) %>%
    select(shoot_biomass_mg, shoot_height, nodule_number) %>%
    # Compute the strain mean
    summarize(across(everything(), function (x) mean(x, na.rm = T))) %>%
    filter(!exp_id == "control") %>%
    ungroup()

write_csv(gst_r, paste0(folder_data, "phenotypes/tradeoff/gst_r.csv"))
write_csv(plants_r, paste0(folder_data, "phenotypes/tradeoff/plants_r.csv"))


# Compute the PCs
gst_pc <- gst_r[,-1] %>% prcomp(center = TRUE, scale = TRUE)
gst_pcs <- as_tibble(summary(gst_pc)$x) %>% mutate(exp_id = gst_r$exp_id) %>% select(exp_id, everything())
gst_pc_imp <- as_tibble(as.list(summary(gst_pc)$importance[2, ])) # proportion of variance


temp <- plants_r %>%
    select(-genome_id, -gradient, -population) %>%
    group_by(exp_plant, exp_nitrogen) %>%
    nest() %>%
    mutate(
        data_rmna = map(data, ~ select(.x, where(~ !all(is.na(.))))),
        exp_ids = map(data_rmna, ~select(.x, exp_id)),
        plants_pc = map(data_rmna, ~ prcomp(select(.x, -exp_id), center = TRUE, scale. = TRUE)),
        plants_pcs = map(plants_pc, ~ as_tibble(summary(.x)$x)),
        plants_pc_imp = map(plants_pc, ~as_tibble(as.list(summary(.x)$importance[2,])))
    )

plants_pcs <- temp %>%
    unnest(c(plants_pcs, exp_ids)) %>%
    select(exp_plant, exp_nitrogen, exp_id, starts_with("PC"))

plants_pc_imp <- temp %>%
    unnest(plants_pc_imp) %>%
    select(exp_plant, exp_nitrogen, starts_with("PC"))


write_csv(gst_pcs, paste0(folder_data, "phenotypes/tradeoff/gst_pcs.csv"))
write_csv(gst_pc_imp, paste0(folder_data, "phenotypes/tradeoff/gst_pc_imp.csv"))
write_csv(plants_pcs, paste0(folder_data, "phenotypes/tradeoff/plants_pcs.csv"))
write_csv(plants_pc_imp, paste0(folder_data, "phenotypes/tradeoff/plants_pc_imp.csv"))


