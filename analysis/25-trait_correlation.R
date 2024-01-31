#' This script compares the correaltion between the two composite traits: growth and symbiosis

renv::load()
library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

gc_prms <- read_csv(paste0(folder_data, "temp/21-gc_prms.csv"))
plants <- read_csv(paste0(folder_data, "temp/23-plants.csv"))
isolates_mapping
unique(gc_prms$exp_id)
unique(plants$rhizobia)

# 0. Clean up
tr_gc <- gc_prms %>%
    left_join(isolates_mapping) %>%
    mutate(gc_unique_id = seq_len(n())) %>%
    select(exp_id, gc_unique_id, temperature, r, lag, maxOD)

tr_gc_long <- tr_gc %>%
    arrange(exp_id) %>%
    select(-gc_unique_id) %>%
    pivot_longer(-c(exp_id, temperature)) %>%
    mutate(trait = paste0(temperature, "_", name)) %>% 
    select(exp_id, trait, value)

tr_sym <- plants %>%
    left_join(isolates_mapping) %>%
    select(exp_id, sym_unique_id = id, contains("mg"), contains("nodule"), contains("px")) 

tr_sym_long <- tr_sym %>%
    arrange(exp_id) %>%
    select(-sym_unique_id) %>%
    pivot_longer(-exp_id, names_to = "trait")



# 1. Find the PCs for growth and symbiosis traits, respectively
# 1.1 Take average across techinical replicates
tr_gc_wide <- tr_gc %>%
    group_by(exp_id, temperature) %>%
    summarize(r = mean(r, na.rm = T), lag = mean(r,na.rm = T), maxOD = mean(maxOD,na.rm = T)) %>%
    pivot_wider(names_from = temperature, values_from = c(r, lag, maxOD)) %>%
    ungroup()
tr_sym_wide <- tr_sym %>%
    group_by(exp_id) %>%
    summarize(nodule_number = mean(nodule_number, na.rm = T), dry_weight_mg = mean(dry_weight_mg, na.rm = T)) %>%
    filter(exp_id != "filter")

# 1.2 Calculate PCs
tr_gc_wide[is.na(tr_gc_wide)] <- 0 # This is not correct
pca_gc <- prcomp(tr_gc_wide[,-1])
tr_gc_pca <- tr_gc_wide[,1] %>% bind_cols(pca_gc$x) %>%
    rename_with(function(x) {paste0(x, "_gc")}, starts_with("PC"))

pca_sym <- prcomp(tr_sym_wide[,-1])
tr_sym_pca <- tr_sym_wide[,1] %>% bind_cols(pca_sym$x) %>%
    rename_with(function(x) {paste0(x, "_sym")}, starts_with("PC"))

# 1.3 Bind symbionsis and growth traits
tr_pca <- tr_sym_pca %>% 
    left_join(tr_gc_pca) %>%
    drop_na()




# 2. Permute the trait pair
# 2.1 Find the strain intersection between the two experiments
list_strains <- intersect(tr_gc_long$exp_id, tr_sym_long$exp_id)
tr_long <- bind_rows(tr_gc_long, tr_sym_long)
tr_long <- tr_long %>%
    filter(exp_id %in% list_strains)

# 2.2 permute the data 
# Create a list of trait combinations 
tb_trait <- crossing(trait1 = unique(tr_gc_long$trait), trait2 = unique(tr_sym_long$trait))
tb_trait$boots <- NA

n_draws <- 100 # Nnmber of bootstrap
n_reps <- 100 # Number of trait pairs to draw in one bootstrap

for (k in seq_len(nrow(tb_trait))) {

    trait1 <- tb_trait$trait1[k]
    trait2 <- tb_trait$trait2[k]

    tb_cors <- tibble(bootstrap = seq_len(n_draws), cor_within = NA, cor_across = NA)

    print(k)
    for (i in seq_len(n_draws)) {
        set.seed(i)
        tr1 <- tr_long %>%
            filter(trait == trait1) %>%
            rename(trait1 = value) %>%
            select(-trait)

        tr2 <- tr_long %>%
            filter(trait == trait2) %>%
            rename(trait2 = value) %>%
            select(-trait)

        # Create a table of all within-strain pairs of traits
        trs <- full_join(tr1, tr2, relationship = "many-to-many", by = join_by(exp_id))

        # Draw random pairs within strains
        trs_within <- trs[sample(seq_len(nrow(trs)), n_reps, replace = T),c("trait1", "trait2")]
        trs_within <- trs_within %>% drop_na()

        # Draw random pairs across strains
        trs_across <- tibble(
            trait1 = tr1$trait1[sample(seq_len(nrow(tr1)), n_reps, replace = T)],
            trait2 = tr2$trait2[sample(seq_len(nrow(tr2)), n_reps, replace = T)]
        )
        trs_across <- trs_across %>% drop_na()

        # Compare correlation
        tb_cors$cor_within[i] <- cor(trs_within$trait1, trs_within$trait2, method = "spearman")
        tb_cors$cor_across[i] <- cor(trs_across$trait1, trs_across$trait2, method = "spearman")
    }

    tb_trait$boots[k] <- list(tb_cors)

}

# 2.3 Bind the traits
tb_traits <- tb_trait %>% unnest(boots)

tb_traits <- tb_traits %>%
    group_by(trait1, trait2) %>%
    mutate(dif = case_when(
        cor_within > cor_across ~ "greater",
        cor_within < cor_across ~ "less"
    ))


write_csv(tb_traits, paste0(folder_data, "temp/24-tb_traits.csv"))

# 
tb_traits_dif <- tb_traits %>%
    mutate(dif = factor(dif, c("greater", "less"))) %>%
    group_by(trait1, trait2, dif, .drop = F) %>%
    count() 

write_csv(tb_traits_dif, paste0(folder_data, "temp/24-tb_traits_dif.csv"))


tb_traits_dif_sig <- tb_traits_dif %>% filter(n >= 95) 
write_csv(tb_traits_dif_sig, paste0(folder_data, "temp/24-tb_traits_dif_sig.csv"))    



# 2.4 Save the drawn pair data of one trait pair
generate_reps <- function (trait1, trait2) {

tb_cors <- tibble(bootstrap = seq_len(n_draws), cor_within = NA, cor_across = NA)
tb_rep <- rep(list(NA), n_reps)

for (i in seq_len(n_draws)) {
    set.seed(i)
    tr1 <- tr_long %>%
        filter(trait == trait1) %>%
        rename(trait1 = value) %>%
        select(-trait)

    tr2 <- tr_long %>%
        filter(trait == trait2) %>%
        rename(trait2 = value) %>%
        select(-trait)

    # Create a table of all within-strain pairs of traits
    trs <- full_join(tr1, tr2, relationship = "many-to-many", by = join_by(exp_id))

    # Draw random pairs within strains
    trs_within <- trs[sample(seq_len(nrow(trs)), n_reps, replace = T),c("trait1", "trait2")]
    trs_within <- trs_within %>% drop_na()

    # Draw random pairs across strains
    trs_across <- tibble(
        trait1 = tr1$trait1[sample(seq_len(nrow(tr1)), n_reps, replace = T)],
        trait2 = tr2$trait2[sample(seq_len(nrow(tr2)), n_reps, replace = T)]
    )
    trs_across <- trs_across %>% drop_na()

    # Compare correlation
    tb_rep[[i]] <- bind_rows(
        tibble(pair = "within", trait1 = trs_within$trait1, trait2 = trs_within$trait2),
        tibble(pair = "across", trait1 = trs_across$trait1, trait2 = trs_across$trait2),
    ) %>% mutate(bootstrap = i)
    tb_cors$cor_within[i] <- cor(trs_within$trait1, trs_within$trait2)
    tb_cors$cor_across[i] <- cor(trs_across$trait1, trs_across$trait2)
}

    tb_reps <- bind_rows(tb_rep) %>%
        rename(trait1_value = trait1, trait2_value = trait2) %>%
        mutate(trait1 = trait1, trait2 = trait2)
    return(tb_reps)
}

# Give two examples
tb_reps_exp1 <- generate_reps(trait1 = "30c_r", trait = "nodule_number") # THis is significant
tb_reps_exp2 <- generate_reps(trait1 = "35c_r", trait = "nodule_number") # THis is not significant

write_csv(tb_reps_exp1, paste0(folder_data, "temp/24-tb_reps_exp1.csv"))
write_csv(tb_reps_exp2, paste0(folder_data, "temp/24-tb_reps_exp2.csv"))


# 3. Calcualte the mean values 
tr_gc_wide <- tr_gc %>%
    group_by(exp_id, temperature) %>%
    summarize(r_se = sd(r, na.rm = T), lag_se = sd(r,na.rm = T), maxOD_se = sd(maxOD,na.rm = T),
                r = mean(r, na.rm = T), lag = mean(r,na.rm = T), maxOD = mean(maxOD,na.rm = T)) %>%
    pivot_wider(names_from = temperature, values_from = c(r, lag, maxOD, r_se, lag_se, maxOD_se)) %>%
    ungroup()
tr_sym_wide <- tr_sym %>%
    group_by(exp_id) %>%
    summarize(nodule_number_se = sd(nodule_number, na.rm = T), dry_weight_mg_se = sd(dry_weight_mg, na.rm = T),
                nodule_number = mean(nodule_number, na.rm = T), dry_weight_mg = mean(dry_weight_mg, na.rm = T)) %>%
    filter(exp_id != "filter")

tr_traits_mean <- tr_gc_wide %>% left_join(tr_sym_wide)
write_csv(tr_traits_mean, paste0(folder_data, "temp/24-tb_traits_mean.csv"))

