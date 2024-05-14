library(tidyverse)
library(syntenet)
Sys.setenv(PATH = paste("/opt/homebrew/bin", Sys.getenv("PATH"), sep = ":")) # Add homebrew command to the PATH

create_species_id_table <- function(species_names) {

    # Replace space with 'X'
    list_names <- gsub(" ", "X", species_names)

    # Get the minimum number of characters possible
    for(n in seq(3, 6)) {
        abbrev <- substr(list_names, 1, n)
        count <- table(abbrev)
        if (!any(count > 1)) break
    }

    # Create a data frame of species IDs and names
    abbrev_df <- data.frame(
        species_id = abbrev,
        species_name = species_names
    )

    # If there are duplicated names with 5 characters, append LETTER to the end
    if(n > 5) {
        abbrev <- substr(list_names, 1, 5)
        abbrev_df <- data.frame(
            species_name = species_names,
            species_id = abbrev,
            count = as.numeric(ave(
                as.character(abbrev), abbrev, FUN = seq_along
            ))
        )
        abbrev_df$species_id <- ifelse(
            abbrev_df$count > 1,
            paste0(substr(abbrev_df$species_id, 1, 4), LETTERS[abbrev_df$count]),
            abbrev_df$species_id
        )
        abbrev_df <- abbrev_df[, c("species_id", "species_name")]
    }

    return(abbrev_df)
}

# Protein sequences
data(proteomes)
head(proteomes) # Olucimarinus and OspRCC809
class(proteomes) # list of AAStringSet
sapply(proteomes, length) # 1901 1433

# Annotation
data(annotation)
head(annotation)
class(annotation) # CompressedGRangesList
sapply(annotation, length) # 1903 1433

# # Import data
# fasta_dir <- system.file("extdata", "sequences", package = "syntenet")
# dir(fasta_dir) # "Olucimarinus.fa.gz" "OspRCC809.fa.gz"
# aastringsetlist <- fasta2AAStringSetlist(fasta_dir)
# gff_dir <- system.file("extdata", "annotation", package = "syntenet")
# dir(gff_dir) # "Olucimarinus.gff3.gz" "OspRCC809.gff3.gz"
# grangeslist <- gff2GRangesList(gff_dir)

# Data preprocessing
check_input(proteomes, annotation) # This does the following three checks
names(proteomes); names(annotation) # The names of species match
sapply(proteomes, length) # The number of sequneces in each species is not greater than the number of genes in annotation
names(proteomes$Olucimarinus)[1:10]; annotation$Olucimarinus$gene_id[1:10] # The names of sequence match

# Data processing
pdata <- process_input(proteomes, annotation)

# Check if diamond is installed
diamond_is_installed()

# Run DIAMOND
#data(blast_list)
blast_list <- run_diamond(seq = pdata$seq)
names(blast_list)
head(blast_list$Olucimarinus_Olucimarinus)

# Infer synteny network
net <- infer_syntenet(blast_list, pdata$annotation)
class(net)
head(net)

# Get a 2-column data frame of species IDs and names
id_table <- create_species_id_table(names(proteomes))
id_table

# Cluster network
data(network)
head(network)
clusters <- cluster_network(network)
head(clusters)

# Phylogenomic profiling
profiles <- phylogenomic_profile(clusters)
head(profiles)

plot_profiles(profiles)


# Microsynteny-based phylogeny reconstruction
bt_mat <- binarize_and_transpose(profiles)
# Looking at the first 5 rows and 5 columns of the matrix
bt_mat[1:5, 1:5]


if (F) {

#cmd = "/opt/homebrew/bin/diamond"
cmd = "diamond"
args = "help"

tryCatch(
    system2(cmd, args = args, stdout = FALSE, stderr = FALSE),
    error = function(e) return(FALSE),
    warning = function(w) return(FALSE)
)
}
