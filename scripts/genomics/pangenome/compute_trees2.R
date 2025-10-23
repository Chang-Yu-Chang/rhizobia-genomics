#' Compute trees based on gene content variation

library(tidyverse)
library(ape)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

compute_gpa_tree <- function () {
    tt <- read_gpas()
    d <- dist(t(tt$gpa[,-1]))
    tr <- hclust(d, method = "ward.D2") %>% as.phylo()
    return(tr)
}
compute_spa_tree <- function () {
    tt <- read_gpas()
    d <- dist(t(tt$spa[,-1]))
    tr <- hclust(d, method = "ward.D2") %>% as.phylo()
    return(tr)
}

tr_gpa <- compute_gpa_tree()
tr_spa <- compute_spa_tree()

write.tree(tr_gpa, paste0(folder_genomics, "pangenome/trees/tr_gpa.tree"))
write.tree(tr_spa, paste0(folder_genomics, "pangenome/trees/tr_spa.tree"))
