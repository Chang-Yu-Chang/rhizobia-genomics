#' This script matches the plant phenotypes into one tabular csv

library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 1. Merge the files
treatments <- read_csv(paste0(folder_data, "raw/rhizobia/04-manual_phenotyping/treatments_assigned.csv"), show_col_types = F) %>%
    rename(dry_weight = `DryWeight (mg)`) %>%
    clean_names()

features <- read_csv(paste0(folder_data, "raw/rhizobia/05-root_architecture/features.csv"), show_col_types = F) %>%
    clean_names() %>%
    mutate(id = str_replace(file_name, ".png", "") %>% as.numeric())
treatments <- treatments %>% left_join(features)

traits <- c("dry_weight", "nodule_number", "number_of_root_tips", "number_of_branch_points",
            "total_root_length_px", "branching_frequency_per_px", "network_area_px2",
            "average_diameter_px", "median_diameter_px", "maximum_diameter_px",
            "perimeter_px", "volume_px3", "surface_area_px2")

treatments_long <- treatments %>%
    select(-contains("range_")) %>%
    pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value")


write_csv(treatments_long, paste0(folder_data, "temp/03-treatments_long.csv"))
write_csv(treatments, paste0(folder_data, "temp/03-treatments.csv"))


# 2. Check assumptions
treatments_scaled <- treatments %>%
    # # Excluding the strains that do not nodulate
    # filter(!rhizobia %in% c("H2M3R1", "L4M2R2")) %>%
    # scale the traits
    mutate_at(c(9, 10, 13:48), ~ c(scale(.)))

treatments_scaled_long <- treatments_scaled %>%
    select(-contains("range_")) %>%
    pivot_longer(cols = all_of(traits), names_to = "trait", values_to = "value")

write_csv(treatments_scaled, paste0(folder_data, "temp/03-treatments_scaled.csv"))
write_csv(treatments_scaled_long, paste0(folder_data, "temp/03-treatments_scaled_long.csv"))





