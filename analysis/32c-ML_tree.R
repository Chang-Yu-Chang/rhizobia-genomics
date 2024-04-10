#' This script make the distance matrix for NJ trees

renv::load()
library(tidyverse)
library(janitor)
library(ape)
library(phangorn)
library(tidytree)
source(here::here("analysis/00-metadata.R"))

#tr <- read.tree(paste0(folder_data, "genomics/mltree/aarA/aln.treefile"))
tr <- read.tree(paste0(folder_data, "genomics/mltree/core_b/aln.treefile"))
plot(tr)
class(tr)
tr$edge.length
tr$node.label
tr$tip.label


if (FALSE) {
dists <- read_csv(paste0(folder_data, 'temp/31-dists.csv'))
isolates <- read_csv(paste0(folder_data, "temp/00-isolates.csv"))
    # Compute gene content ML tree
    gpa <- read_csv(paste0(folder_data, "temp/13-gpa.csv"))

    set.seed(1)
    tb <- tibble(x = rbinom(10, 1, 0.5), y = rbinom(10, 1, 1), z = rbinom(10, 1, 0.3))

    dist.ml(tb)

    ace(tb, type = "discrete", method = "ML")

    fdir <- system.file("extdata/trees", package = "phangorn")
    primates <- read.phyDat(file.path(fdir, "primates.dna"), format = "interleaved")
    class(primates)
    dm <- dist.ml(primates); class(dm)
    treeUPGMA <- upgma(dm); class(treeUPGMA) # unweighted pair group method with arithmetic mean
    treeNJ <- NJ(dm)
    plot(treeUPGMA, main="UPGMA")
    plot(treeNJ, "unrooted", main="NJ")

    fun <- function(x) upgma(dist.ml(x))
    bs_upgma <- bootstrap.phyDat(primates, fun); class(bs_upgma)
    plotBS(treeUPGMA, bs_upgma, main="UPGMA")

    #
    parsimony(treeUPGMA, primates)
    parsimony(treeNJ, primates)

    #
    treeRatchet  <- pratchet(primates, trace = 0, k = 10, minit = 1000)
    parsimony(treeRatchet, primates)

    # assign edge length (number of substitutions)
    treeRatchet <- acctran(treeRatchet, primates)

    treeRatchet <- di2multi(treeRatchet)

    class(treeRatchet)
    plotBS(midpoint(treeRatchet), type = "phylogram")
    add.scale.bar()

    treeRA <- random.addition(primates)
    treeSPR  <- optim.parsimony(treeRA, primates)


    #
    mt <- modelTest(primates)

}



