#

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
g <- igraph::read_graph(paste0(folder_data, "genomics/pangenome/all/final_graph.gml"), format = "gml")

# Save nodules and edges
write_graph(g, paste0(folder_data, "phylogenomics_analysis/synteny/edgelist.txt"), format = "edgelist")
nodes <- tibble(id = 1:length(V(g)), gene = V(g)$name)
write_delim(nodes, delim = "\t", paste0(folder_data, "phylogenomics_analysis/synteny/nodes.txt"))

#
# #
# tt <- read_gpas()
# nodes <- read_delim(delim = "\t", paste0(folder_data, "phylogenomics_analysis/synteny/nodes.txt")) %>%
#     filter(str_detect(gene, "nif|nod|fix"))
# edges <- read_table(paste0(folder_data, "phylogenomics_analysis/synteny/edgelist.txt"), col_names = F) %>%
#     rename(from = X1, to = X2) %>%
#     filter(from %in% nodes$id, to %in% nodes$id)
#
# range(nodes$id)
# range(edges$from)
# range(edges$to)
#
# gg <- tbl_graph(
#     nodes = nodes,
#     edges = edges
# )
# THE GRAPH TO TOO LARGE TO MAKE IN GGRAPH



# g_test <- make_ring(3)
# write_graph(g_test, "test.gml", format = "gml")
#
# g <- read_graph("test.gml", format = "gml")


# tbl_graph(
#     nodes = tibble(id = 1:10),
#     edges = tibble(from = c(1,3,5), to = c(4, 9, 3))
# )
