# This script stores the metadata shared by all scripts

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(RColorBrewer)
})

# This main folder depends on your home directory and user name
folder_script <- "~/Desktop/lab/local-adaptation/analysis/" # Enter the directory of analysis scripts
folder_data <- "~/Dropbox/lab/local-adaptation/data/" # Enter the directory of data
folder_genomes <- paste0(folder_data, "genomics/genomes/")

# Table for genomic workflow
genomes_mapping <- tibble(
    batch_name = c(rep("Chang_Q5C_results", 10), "Chang_Q5C_results_repeated", rep("Chang_Q5C_results", 6), "Chang_Q5C_results_repeated", "Chang_Q5C_results",
                   rep("Chang_W8S_results", 18),
                   rep("ncbi", 4)),
    sample_id = c(paste0("Chang_Q5C_", 1:19), paste0("Chang_W8S_", 1:18), "usda1106", "em1021", "em1022", "wsm419"),
    accession = c(rep(NA, 37), c("GCF_002197065.1", "GCF_000006965.1", "GCF_013315775.1", "GCF_000017145.1")),
    genome_name = c(paste0("Chang_Q5C_", 1:19), paste0("Chang_W8S_", 1:18), "usda1106", "em1021", "em1022", "wsm419"),
    genome_id = factor(c(paste0("g", 1:37), "usda1106", "em1021", "em1022", "wsm419")),
)
write_csv(genomes_mapping, paste0(folder_data, "temp/00-genomes_mapping.csv"))

# Table for the 32 strains with growth curve and/or plant inoculation data
isolates_mapping <- tibble(
    exp_id = c("H2M3R1", "H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1",
               "L2M4R1", "L3M5R1", "L4M2R2", "L4M3R3", "L4M4R1", "L4M7R1", "L3M1R1", "src-2 ",
               "fp1-2", "fp1-3", "fp2-2", "crp1-2", "crp1-3", "crp2-2", "gp1-1", "gp1-2",
               "gp1-3", "bg-2", "bg-3", "pms-1", "pms-2", "pms-3", "ppf-1", "40th-1"),
    rhizobia_site = c(rep("high-elevation", 6), rep("low-elevation", 9),
                      rep("suburban", 10), rep("urban", 7)),
    sample_id = c(paste0("Chang_Q5C_", c(2:6, 8:11, 13, 15:17, 19)), paste0("Chang_W8S_", c(1:18))),
    genome_name = c(paste0("Chang_Q5C_", c(2:6, 8:11, 13, 15:17, 19)), paste0("Chang_W8S_", c(1:18))),
    genome_id = paste0("g", c(2:6, 8:11, 13, 15:17, 19:37))
)

write_csv(isolates_mapping, paste0(folder_data, "temp/00-isolates_mapping.csv"))

# Table for both
isolates <- full_join(genomes_mapping, isolates_mapping) %>%
    # Fill in for ncbi strains
    mutate(rhizobia_site = ifelse(batch_name == "ncbi", "ncbi", rhizobia_site),
           exp_id = ifelse(batch_name == "ncbi", "ncbi", exp_id)) %>%
    select(genome_name, genome_id, exp_id, rhizobia_site)

write_csv(isolates, paste0(folder_data, "temp/00-isolates.csv"))

# Trait vectors
traits <- c("dry_weight_mg", "nodule_number", "root_weight_mg", "nodule_weight_mg",
            "number_of_root_tips", "number_of_branch_points",
            "total_root_length_px", "branching_frequency_per_px", "network_area_px2",
            "average_diameter_px", "median_diameter_px", "maximum_diameter_px",
            "perimeter_px", "volume_px3", "surface_area_px2")
traits2 <- c(traits,
             paste0(rep(c("root_length_diameter_range_", "projected_area_diameter_range_", "surface_area_diameter_range_", "volume_diameter_range_"), each = 6),
                    rep(1:6, 4), rep(c("_px", "_px2", "_px2", "_px3"), each = 6)))


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


# Plotting functions
paint_white_background <- function () {
    theme(plot.background = element_rect(color = NA, fill = "white"))
}

theme_facets <- theme(
    panel.grid.major = element_line(color = "grey95"),
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.border = element_rect(fill = NA, color = "black")
)


# Colors
rhizobia_site_colors <- brewer.pal(n = 6, name = "Paired")[c(1,2,5,6)] %>% setNames(c("high-elevation", "low-elevation", "suburban", "urban"))






# list_strains <- c("H1M1R1", "H2M3R1", "H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H3M4R2", "H4M5R1", "L1M2R2", "L2M2R1", "L2M4R1", "L3M4R1", "L3M5R1", "L3M6R2", "L4M2R2", "L4M3R3", "L4M4R1", "L4M4R3", "L4M7R1", "blank")
# rhizobia_strains <- c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2")
# rhizobia_alphas <- setNames(c(.5,.7,.9, .5,.7,.9, .5), c("H2M3R1", "H3M1R1", "H4M5R1", "L2M2R1", "L3M5R1", "L4M2R2", NA))
# rhizobia_site_colors <- c(`high-elevation` = "#0C6291", `mid-elevation` = "#CBD4C2", `low-elevation` = "#BF4342", `urban` = "deeppink", `suburban` = "yellowgreen")
# rhizobia_site_colors <- c(`high-elevation` = "#1f78b4", `low-elevation` = "#ff7f00", `urban` = "#33a02c", `suburban` = "#6a3d9a")
