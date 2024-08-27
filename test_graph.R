
library(igraph)
renv::load()
library(tidyverse)
library(gggenomes)
library(cowplot)
source(here::here("metadata.R"))

gg <- read_graph(paste0(folder_data, "genomics/pangenome/isolates/final_graph.gml"), format = "gml")
