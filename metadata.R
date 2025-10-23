#' This script stores the metadata shared by all R scripts

# This main folder depends on your home directory and user name
folder_data <- "~/Dropbox/lab/rhizobia-genomics/data/" # Enter the directory of data
folder_genomics <- paste0(folder_data, "genomics/")
folder_phenotypes <- paste0(folder_data, "phenotypes/")

# Palettes for plotting
species_colors <- c(`Ensifer adhaerens` = "#B67C66", `Ensifer sp.` = "grey20", `Sinorhizobium medicae` = "#4DCCBD", `Sinorhizobium meliloti` = "#B07BAC")
species_shapes <- c(`Ensifer adhaerens` = 15, `Ensifer sp.` = 16, `Sinorhizobium medicae` = 22, `Sinorhizobium meliloti` = 21)

# utils
read_gpas <- function () {
    folder_gpas <- paste0(folder_genomics, "pangenome/gene_content/")
    isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
    gpa <- read_csv(paste0(folder_gpas, "gpa.csv"))
    gpar <- read_csv(paste0(folder_gpas, "gpar.csv"))
    list_sccg <- read_csv(paste0(folder_gpas, "list_sccg.csv"), col_names = "gene")
    sml <- read_csv(paste0(folder_gpas, "sml.csv"))
    spa <- read_csv(paste0(folder_gpas, "spa.csv"))
    gpatl <- read_csv(paste0(folder_gpas, "gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))
    gd <- read_csv(paste0(folder_gpas, "gd.csv"))
    # gene_order <- read_csv(paste0(folder_gpas, "gene_order.csv"))
    # gpacl <- read_csv(paste0(folder_gpas, "gpacl.csv")) %>%
    #     mutate(genome_id = factor(genome_id, rev(isolates$genome_id)))
    gcn <- read_csv(paste0(folder_gpas, "gcn.csv"))
    cleaned_gene_names <- read_csv(paste0(folder_gpas, "cleaned_gene_names.csv"))
    return(list(gpa = gpa, gpar = gpar, list_sccg = list_sccg, sml = sml, spa = spa, gpatl = gpatl,
                #gene_order = gene_order, gpacl = gpacl,
                gd = gd, gcn = gcn, cleaned_gene_names = cleaned_gene_names))
}
# turn_p_to_asteriks <- function (p_value) {
#     if (p_value < 0.001) {
#         asterisks <- "***"
#     } else if (p_value < 0.01) {
#         asterisks <- "**"
#     } else if (p_value < 0.05) {
#         asterisks <- "*"
#     } else {
#         asterisks <- "n.s."
#     }
#     return(asterisks)
# }
#
# clean_p_lab <- function (p_value) {
#         p_value = round(p_value, 3)
#     if (p_value < 0.001) {
#         cp <- "p<0.001***"
#     } else if (p_value < 0.01) {
#         cp <- paste0("p=", p_value, " **")
#     } else if (p_value < 0.05) {
#         cp <- paste0("p=", p_value, " *")
#     } else {
#         cp <- paste0("p=", p_value, " n.s.")
#     }
#     return(cp)
# }
#
# # Flextable
# edit_p <- function (pv) {
#     if (pv < 0.001) {
#         return("<0.001")
#     } else {
#         return(as.character(round(pv, 3)))
#     }
# }
# detect_sig <- function (pv) {
#     if (pv > 0.05) {
#         return("n.s.")
#     } else if (pv > 0.01) {
#         return("*")
#     } else if (pv > 0.001) {
#         return("**")
#     } else if (pv < 0.001) {
#         return("***")
#     }
# }
