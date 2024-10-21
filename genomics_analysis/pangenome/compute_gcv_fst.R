#' This script plots the heatmap of gene content

library(tidyverse)
library(janitor)
library(ape) # for reading multiple fasta
library(poppr) # for reading snps into a genind object
library(hierfstat) # for computing Fst
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
make_genind_from_gcv <- function (gpa, gene_name) {
    gpa_i <- t(gpa[which(gpa$gene == gene_name),-1])
    pop <- isolates$population[match(rownames(gpa_i), isolates$genome_id)]
    df2genind(gpa_i, pop = pop, ploidy = 1)
}


for (set_name in c("elev_med", "urbn_mel")) {
    #set_name = "elev_med"
    #set_name = "urbn_mel"

    tt <- read_gpas(set_name)
    gcv_fst <- list()
    for (gene in tt$gpa$gene) {
        #gene = "xerC_11~~~xerC_5"
        genind_object <- make_genind_from_gcv(tt$gpa, gene)
        if (length(genind_object@loc.fac) == 1) {cat(gene, "is a core gene \n"); next}
        ww <- wc(genind_object, diploid = F)
        gcv_fst[[gene]] <- tibble(fst = ww$FST)
        cat("\nProcessed:", gene)
    }
    per_acce_fst <- bind_rows(gcv_fst, .id = "gene") # per acce gene F_st results
    write_csv(per_acce_fst, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/per_acce_fst.csv"))
    # top_gene_or <- slice_top_genes(tidied_fisher)
    # write_csv(top_gene_or, paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv"))

}
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
    slice_top_genes <- function (tidied_fisher, prop = 0.01) {
        # Genes with top ORs
        tidied_fisher %>%
            ungroup() %>%
            arrange(desc(or)) %>%
            slice_max(or, prop = prop) %>% # The top 1%
            filter(!str_detect(gene, "group")) # Remove gene annotated as "group_XX"

    }

}




