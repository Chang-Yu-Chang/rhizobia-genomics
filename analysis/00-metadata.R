# This script stores the metadata shared by all scripts

library(tidyverse)

# This main folder depends on your home directory and user name
folder_script <- "~/Desktop/lab/local-adaptation/analysis/" # Enter the directory of analysis scripts
folder_data <- "~/Dropbox/lab/local-adaptation/data/" # Enter the directory of data

# For genomics
genomes <- tibble(genome_name = c(paste0("Chang_Q5C_", c(2:6, 8:11, 13, 15:17, 19)), "em1021", "em1022", "wsm419"),
                  genome_id = factor(paste0("g", 1:17)))


# FOr experiemnt
list_strains <- c("H1M1R1", "H2M3R1", "H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H3M4R2", "H4M5R1", "L1M2R2", "L2M2R1", "L2M4R1", "L3M4R1", "L3M5R1", "L3M6R2", "L4M2R2", "L4M3R3", "L4M4R1", "L4M4R3", "L4M7R1", "blank")
rhizobia_strains <- c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")
rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2", NA))
rhizobia_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
plant_site_colors <- c(H = "#0C6291", S = "#CBD4C2", L = "#BF4342")
traits <- c("dry_weight_mg", "nodule_number", "root_weight_mg", "nodule_weight_mg",
            "number_of_root_tips", "number_of_branch_points",
            "total_root_length_px", "branching_frequency_per_px", "network_area_px2",
            "average_diameter_px", "median_diameter_px", "maximum_diameter_px",
            "perimeter_px", "volume_px3", "surface_area_px2")
traits2 <- c(traits,
             paste0(rep(c("root_length_diameter_range_", "projected_area_diameter_range_", "surface_area_diameter_range_", "volume_diameter_range_"), each = 6),
                    rep(1:6, 4), rep(c("_px", "_px2", "_px2", "_px3"), each = 6)))


rhizobia <- tibble(
    strain = list_strains[-20],
    rhizobia_site = str_sub(list_strains[-20], 1, 1)
)


paint_white_background <- function () {
    theme(plot.background = element_rect(color = NA, fill = "white"))
}

trait_axis_names <- c(
    "dry_weight_mg" = "shoot biomass (mg)",
    "nodule_number" = "number of nodules",
    "root_weight_mg" = "root biomass (mg)",
    "nodule_weight_mg" = "nodule biomass (mg)",
    "number_of_root_tips" = "number of root tips",
    "number_of_branch_points" = "number of branch points",
    "total_root_length_px" = "root length (px)",
    "branching_frequency_per_px" = "branching frequencing (1/px)",
    "network_area_px2" = "root area (px^2)",
    "average_diameter_px" = "average diameter (px)",
    "median_diameter_px" = "median diameter (px)",
    "maximum_diameter_px" = "maximum diameter (px)",
    "perimeter_px" = "perimeter (px)",
    "volume_px3" = "volume (px^3)",
    "surface_area_px2" = "surface area (px^2)"
)
