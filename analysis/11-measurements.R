#' This script analyses the hand measured phenotypes: dry weight, nodule count, root weight, nodule weight

library(tidyverse)
library(broom)
library(janitor)
# library(lattice)
# library(lsmeans)
# library(glmmADMB)
source(here::here("analysis/00-metadata.R"))

# Manual measurements
treatments <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(dry_weight = `DryWeight (mg)`) %>%
    clean_names()
nrow(treatments) # 167 plants
treatments <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    clean_names()

# Computer vision measurements
features <- read_csv(paste0(folder_data, "raw/rhizobia/05-root_architecture/features.csv"), show_col_types = F) %>%
    clean_names() %>%
    mutate(id = str_replace(file_name, ".png", "") %>% as.numeric())
treatments <- treatments %>% left_join(features)

treatments_long <- treatments %>%
    select(-contains("range_")) %>%
    pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value")


write_csv(treatments, paste0(folder_data, "temp/11-treatments.csv"))
write_csv(treatments_long, paste0(folder_data, "temp/11-treatments_long.csv"))


# Check assumptions ----
treatments_scaled <- treatments %>%
    # # Excluding the strains that do not nodulate
    # filter(!rhizobia %in% c("H2M3R1", "L4M2R2")) %>%
    # scale the traits
    mutate_at(c(9:11, 14:49), ~ c(scale(.)))

treatments_scaled_long <- treatments_scaled %>%
    select(-contains("range_")) %>%
    pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value")

write_csv(treatments_scaled, paste0(folder_data, "temp/11-treatments_scaled.csv"))
write_csv(treatments_scaled_long, paste0(folder_data, "temp/11-treatments_scaled_long.csv"))


