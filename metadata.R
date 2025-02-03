#' This script stores the metadata shared by all scripts

# This main folder depends on your home directory and user name
folder_script <- "~/Desktop/lab/rhizobia-genomics/analysis/" # Enter the directory of analysis scripts
folder_data <- "~/Dropbox/lab/rhizobia-genomics/data/" # Enter the directory of data
folder_genomics <- paste0(folder_data, "genomics/")
folder_phenotypes <- paste0(folder_data, "phenotypes/")

# Palettes for plotting
gradient_colors <- c("PA" = "gold2", "VA" = "olivedrab")
population_colors <- c(`high elevation` = "#0C6291", `low elevation` = "#BF4342", `suburban` = "#0cc45f", `urban` = "#a642bf", control = "grey")
pops_colors <- c(
    `high-high` = "#0C6291",
    `high-low` = "grey50",
    `low-low` = "#BF4342",
    `suburban-suburban` = "#0cc45f",
    `suburban-urban` = "grey50",
    `urban-urban` = "#a642bf"
)
species_colors <- c(adhaerens = "grey", canadensis = "grey", medicae = "steelblue", meliloti = "maroon")
species_shapes <- c(meliloti = 21, medicae = 22, adhaerens = 15, canadensis = 16)
plant_colors <- c(sativa = "#62216d", lupulina = "#fde900")

# Traits
clean_trait_names <- function (x) str_split(x, pattern = " ")[[1]] %>% str_sub(1,1) %>% paste(collapse = "") %>% toupper() %>% str_pad(width = 4, side = "right", pad = " ")
traits <- tibble(
    trait_type = c("nodule", "nodule", "nodule", "nodule", "root", "shoot", "nodule", "leaf", "leaf", "leaf", "nodule", "shoot", "root", "root", "root"),
    trait = c("nodule_number", "nodule_shape", "nodule_size", "nodule_color", "root_biomass_mg", "shoot_biomass_mg", "lateral_root_nodule_number", "leaf_color", "leaf_number", "longest_petiole_length", "primary_root_nodule_number", "shoot_height", "lateral_root_number", "longest_lateral_root_length", "primary_root_length")
) %>%
    mutate(trait_pre = trait %>% str_remove("_mg") %>% str_replace_all("_", " ")) %>%
    #mutate(trait_abr = map_chr(trait_pre, clean_trait_names)) %>%
    mutate(trait_type = factor(trait_type, c("shoot", "nodule", "leaf", "root"))) %>%
    arrange(trait_type) %>%
    mutate(trait_pre = case_when(
        trait_pre == "root biomass" ~ "root biomass (mg)",
        trait_pre == "shoot biomass" ~ "shoot biomass (mg)",
        trait_pre == "longest petiole length" ~ "longest petiole length (cm)",
        trait_pre == "shoot height" ~ "shoot height (cm)",
        trait_pre == "primary root length" ~ "primary root length (cm)",
        trait_pre == "longest lateral root length" ~ "longest lateral root length (cm)",
        #trait_pre == "lateral root number" ~ "lateral\nroot\nnumber",
        T ~ trait_pre
    )) %>%
    mutate(trait_pre2 = factor(str_replace_all(trait_pre, " ", "\n"))) %>%
    mutate(trait_pre = factor(trait_pre, trait_pre))

traits <- bind_rows(traits, tibble(
    trait_type = "growth", trait = c("r", "lag", "maxOD"), trait_pre = c("growth rate (1/hr)", "lag time (hr)", "yield [OD]")
))



# utils
read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gpar <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpar.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))
    gcn <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gcn.csv"))
    cleaned_gene_names <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/cleaned_gene_names.csv"))
    return(list(gpa = gpa, gpar = gpar, list_sccg = list_sccg, sml = sml, spa = spa, gpatl = gpatl, gene_order = gene_order, gd = gd, gpacl = gpacl, gcn = gcn, cleaned_gene_names = cleaned_gene_names))
}
read_fsts <- function (set_name) {
    per_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_gene_fst.csv"))
    per_locus_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/per_locus_fst.csv"))
    gene_lengths <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name,"/gene_lengths.csv"))
    return(list(per_gene_fst = per_gene_fst, per_locus_fst = per_locus_fst, gene_lengths = gene_lengths))
}
turn_p_to_asteriks <- function (p_value) {
    if (p_value < 0.001) {
        asterisks <- "***"
    } else if (p_value < 0.01) {
        asterisks <- "**"
    } else if (p_value < 0.05) {
        asterisks <- "*"
    } else {
        asterisks <- "n.s."
    }
    return(asterisks)
}

clean_p_lab <- function (p_value) {
        p_value = round(p_value, 3)
    if (p_value < 0.001) {
        cp <- "p<0.001***"
    } else if (p_value < 0.01) {
        cp <- paste0("p=", p_value, " **")
    } else if (p_value < 0.05) {
        cp <- paste0("p=", p_value, " *")
    } else {
        cp <- paste0("p=", p_value, " n.s.")
    }
    return(cp)
}

# Flextable
edit_p <- function (pv) {
    if (pv < 0.001) {
        return("<0.001")
    } else {
        return(as.character(round(pv, 3)))
    }
}
detect_sig <- function (pv) {
    if (pv > 0.05) {
        return("n.s.")
    } else if (pv > 0.01) {
        return("*")
    } else if (pv > 0.001) {
        return("**")
    } else if (pv < 0.001) {
        return("***")
    }
}
