#' This script stores the metadata shared by all scripts

# This main folder depends on your home directory and user name
folder_script <- "~/Desktop/lab/local-adaptation/analysis/" # Enter the directory of analysis scripts
folder_data <- "~/Dropbox/lab/local-adaptation/data/" # Enter the directory of data
folder_genomics <- paste0(folder_data, "genomics/")

# Table for genomics workflow
genomes <- tibble(
    batch_name = c(rep("Chang_Q5C_results", 8), "Chang_Q5C_results_repeated", rep("Chang_Q5C_results", 5), rep("Chang_W8S_results", 18), rep("ncbi", 5)),
    genome_name = c(paste0("Chang_Q5C_", c(2:6,8:11,13,15:17,19)), paste0("Chang_W8S_", 1:18), "usda1106", "em1021", "em1022", "wsm419", "casidaa"),
    genome_id = factor(c(paste0("g", c(2:6,8:11,13,15:17,19:37)), "usda1106", "em1021", "em1022", "wsm419", "casidaa")),
    accession = c(rep(NA, 32), c("GCF_002197065.1", "GCF_000006965.1", "GCF_013315775.1", "GCF_000017145.1", "GCF_000697965.2"))
)

write_csv(genomes, paste0(folder_data, "mapping/genomes.csv"))

# Table for the 32 strains with growth curve and/or plant inoculation data
isolates <- tibble(
    exp_id = c("H2M3R1", "H2M3R2", "H3M1R1", "H3M3R2", "H3M4R1", "H4M5R1", "L1M2R2", "L2M2R1",
               "L2M4R1", "L3M5R1", "L4M2R2", "L4M3R3", "L4M4R1", "L4M7R1", "L3M1R1", "src-2",
               "fp1-2", "fp1-3", "fp2-2", "crp1-2", "crp1-3", "crp2-2", "gp1-1", "gp1-2",
               "gp1-3", "bg-2", "bg-3", "pms-1", "pms-2", "pms-3", "ppf-1", "40th-1"),
    site = c("H2", "H2", "H3", "H3", "H3", "H4", "L1", "L2", "L2", "L3", "L4", "L4", "L4", "L4", "L3", "src", "fp", "fp", "fp", "crp", "crp", "crp", "gp", "gp", "gp", "bg", "bg", "pms", "pms", "pms", "ppf", "40th"),
    site_group = c(rep("high elevation", 6), rep("low elevation", 9), rep("suburban", 10), rep("urban", 7)),
    population = c(rep("VA", 15), rep("PA",17)),
    genome_name = c(paste0("Chang_Q5C_", c(2:6, 8:11, 13, 15:17, 19)), paste0("Chang_W8S_", c(1:18))),
    genome_id = paste0("g", c(2:6, 8:11, 13, 15:17, 19:37))
)

write_csv(isolates, paste0(folder_data, "mapping/isolates.csv"))

# Color
population_colors <- c("PA" = "gold2", "VA" = "olivedrab")
site_group_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", `suburban` = "#0cc45f", `urban` = "#a642bf", control = "grey")
species_colors <- c(adhaerens = "grey", canadensis = "grey", medicae = "steelblue", meliloti = "maroon")
species_shapes <- c(meliloti = 21, medicae = 22, adhaerens = 15, canadensis = 16)
















