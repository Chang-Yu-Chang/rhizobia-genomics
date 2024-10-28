#' This R script is a CLT for distinct the blast results such that each gene has one matched ref seq
#'
#' Usage:
#'
#' arg1  set_name, like elev_med or urbn_mel
#' arg2  ref blast db, like ngr234
#' arg3  input blast results file name
#' arg4  output blast results file name

args <- commandArgs(trailingOnly = TRUE)


set_name <- args[1]
ref <- args[2]
blast_input_file <- args[3]
blast_output_file <- args[4]

bb <- read.table(blast_input_file)
split_bb <- split(bb, bb$V1)
bb <- do.call(rbind, lapply(split_bb, function(x) {
    x <- x[order(x$V3, decreasing = T),]
    return(x[1,])
}))

write.table(bb, blast_output_file, col.names = F, quote = F, sep = "\t")
