#' This script plots the heatmap of gene content

library(tidyverse)
library(purrr)
library(broom)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))

read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpatl = gpatl, gpacl = gpacl, gd = gd, sml = sml, list_sccg = list_sccg, spa = spa))
}
compute_fisher <- function (m) {
    a = m[1,1]
    b = m[1,2]
    c = m[2,1]
    d = m[2,2]
    p = choose(a + c, a) * choose(b + d, b) / choose(sum(m), a+b)

    m = m + 0.5
    a = m[1,1]
    b = m[1,2]
    c = m[2,1]
    d = m[2,2]

    if (a*d > b*c) {
        or = log(a*d/(b*c))
    } else {
        or = log(b*c/(a*d))
    }


    return(list(or = or, p = p))
}
make_long_count <- function (gpa) {
    #' COunt the gene presence absence within population
    gpa %>%
        pivot_longer(-gene, names_to = "genome_id", values_to = "is_present") %>%
        left_join(select(isolates, genome_id, site_group)) %>%
        mutate(is_present = factor(is_present)) %>%
        group_by(gene, site_group, is_present, .drop = F) %>%
        summarize(n = n())
}
map_fisher_by_gene <- function (long_counts) {
    #' Perform fisher exact test for each gene and tidy it
    long_counts %>%
        group_by(gene) %>%
        summarise(
            contingency_table = list(matrix(
                n,
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("population A", "population B"), c("absence", "presence"))
            )),
            fisher = map(contingency_table, compute_fisher)
            #tidied = map(fisher, as_tibble(fisher))
        ) %>%
        rowwise() %>%
        mutate(tidied = list(as_tibble(fisher))) %>%
        select(gene, tidied) %>%
        unnest(cols = tidied)
        #unnest(cols = tidied)
}

set_name = "elev_med"
tt <- read_gpas(set_name)
long_counts <- make_long_count(tt$gpa)
tidied_fisher <- map_fisher_by_gene(long_counts)
write_csv(tidied_fisher, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/tidied_fisher.csv"))



set_name = "urbn_mel"
tt <- read_gpas(set_name)
long_counts <- make_long_count(tt$gpa)
tidied_fisher <- map_fisher_by_gene(long_counts)
write_csv(tidied_fisher, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/tidied_fisher.csv"))


# Check
tt$gpa %>%
    filter(gene == "nylB")
results %>%
    filter(gene == "nylB") %>%
    pull(tidied)
fisher_results$contingency_table[fisher_results$gene == "nylB"]


if (F) {


p=0.1
p^4 * (1-p)^6 + (1-p)^4 * p^6

n_toss <- 10000
x <- rbinom(n_toss, 4, p = p) # the prob of pop A has all presence
y <- rbinom(n_toss, 6, p = p) # the prob of all pop B has no presence at all
sum(x == 4)/n_toss * sum(y == 0)/n_toss
tidied_fisher %>%
    filter(estimate != 0) %>% # Core genes have estimate = 0
    arrange(p.value)



tidied_fisher %>%
    filter(estimate == 0, is.infinite(conf.high)) %>%
    left_join(mutate(tt$list_sccg, is_sccg = T)) # Those not matched have paralogs so they are not in the list of SCCG
}
if (F) {
    m = matrix(c(3,1,0,6), nrow = 2)
    m = matrix(c(0,4,6,0), nrow = 2)
    m = matrix(c(4,0,0,6), nrow = 2)
    m = matrix(c(3,1,0,6), nrow = 2)
    compute_fisher(m)


    fisher.test(m) %>% tidy()
    #(m[1,1] * m[2,2]) / (m[1,2] * m[2,1])
}




