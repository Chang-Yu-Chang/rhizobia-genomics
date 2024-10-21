#' This script finds the GO terms for the top gene clusters

library(tidyverse)
#library(janitor)
#library(ape) # the build in read.gff does not handle the fasta sequence trailing in the gff file
#library(rtracklayer) # for reading gff
#library(UniProt.ws) # for interfacing with uniprot database
#library(UniprotR)
source(here::here("metadata.R"))

isolates <- read_csv(paste0(folder_data, "mapping/isolates.csv"))
read_gpas <- function (set_name) {
    gpa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpa.csv"))
    gene_order <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gene_order.csv"))
    gpar <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpar.csv"))
    gpatl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpatl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gpacl <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gpacl.csv")) %>%
        mutate(genome_id = factor(genome_id, rev(isolates$genome_id)), gene = factor(gene, gene_order$gene))
    gd <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/gd.csv"))
    sml <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/sml.csv"))
    list_sccg <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/list_sccg.csv"), col_names = "gene")
    spa <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/spa.csv"))

    return(list(gpa = gpa, gene_order = gene_order, gpar = gpar, gpatl = gpatl, gpacl = gpacl, gd = gd, sml = sml, list_sccg = list_sccg, spa = spa))
}

#set_name = "urbn_mel"
set_name = "elev_med"
tt <- read_gpas(set_name)

xx <- tt$gpar %>%
    #filter(gene == "aar") %>%
    mutate(across(matches("^g\\d"), ~str_remove_all(.x, "/Users/cychang/Dropbox/lab/rhizobia-genomics/data/genomics/fasta/genomes/"))) %>%
    pivot_longer(cols = -c(gene, non_unique_gene_name, annotation)) %>%
    separate_rows(value, sep = ";") %>%
    drop_na(value)

nrow(xx) # 73661

# Remove genes with no good annotation -> there are ~4647 unique genes have annotation among ~10k pangenes
yy <- xx %>%
    filter(!(str_detect(gene, "group") & is.na(non_unique_gene_name) & annotation == "hypothetical protein")) %>%
    distinct(gene, non_unique_gene_name, annotation, name, .keep_all = T) %>%
    distinct(gene, .keep_all = T)

list_gene <- yy %>%
    filter(!str_detect(gene, "group")) %>%
    select(gene) %>%
    separate_rows("gene", sep = "~~~") %>%
    mutate(gene = str_remove(gene, "_\\d")) %>%
    unique() %>%
    arrange(gene)

list_gene$gene %>% cat() # All 2000 genes with known names
#GetProteinGOInfo("Q05572")



library(UniprotR)
#Read Accessions from csv file , Note : Accessions must be in the first column.
Accessions <-GetAccessionList("https://s3.amazonaws.com/csvpastebin/uploads/9571fa356c67a0c7c95e8431799a051a/Accessions.csv")
#Get Taxanomy Information
TaxaObj <- GetNamesTaxa(Accessions)
#Visualize Chromosomes localization
PlotChromosomeInfo(TaxaObj)
#Visualize protein's gene name as Network
PlotGenesNetwork(TaxaObj)

#Get Gene ontolgy Information
GeneOntologyObj <- GetProteinGOInfo(Accessions)
#tibble(GeneOntologyObj) %>% view
#Plot Biological process information top 10 go terms
PlotGOBiological(GeneOntologyObj, Top = 10)
#Plot molecular function information top 20 go terms
Plot.GOMolecular(GeneOntologyObj, Top = 20)
#Plot subcellualr localization information
Plot.GOSubCellular(GeneOntologyObj)
#Combine Gene ontology plots into one plot
PlotGoInfo(GeneOntologyObj)
#Handy visualization for publications
PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = getwd(), width = 8, height = 5)


if (F) {

availableUniprotSpecies(pattern = "meliloti")

allToKeys(fromName = "UniProtKB_AC-ID") %>% str_subset("Gene")

mapUniProt(
    from = "UniProtKB_AC-ID",
    to = "Gene_Name",
    query = "Q05572",
    verbose = T
)

mapUniProt(
    from = "UniProtKB_AC-ID",
    to = "UniProtKB",
    columns = c("accession", "id"),
    query = list(organism_id = 10090, ids = c('Q7TPG8', 'P63318'))
)


from = "UniProtKB_AC-ID"
to = "UniRef90"
columns = character(0L)
query = "Q05572"
verbose = FALSE
debug = FALSE
paginate = TRUE
pageSize = 500L


library(httr)
.UNIPROT_REST_URL <- "https://rest.uniprot.org/"
.messageDEBUG <- function(url, debug) {
    if (debug)
        message("Hitting: ", url)
    url
}
.getResultsURL <- function(redurl, paginate, debug) {
    if (!paginate) {
        redurl <- gsub(
            "/idmapping/results/", "/idmapping/stream/", redurl, fixed = TRUE
        )
        redurl <- gsub("/results/", "/results/stream/", redurl, fixed = TRUE)
    }
    .messageDEBUG(redurl, debug)
}
.handleResults <- function(results, debug) {
    rdata <- read.delim(text = content(results, encoding = "UTF-8"))
    while (length(headers(results)$link)) {
        nextlink <- headers(results)$link
        results <- GET(
            .messageDEBUG(gsub("<(.*)>.*", "\\1", nextlink), debug),
            accept_json()
        )
        result <- read.delim(text = content(results, encoding = "UTF-8"))
        rdata <- do.call(rbind.data.frame, list(rdata, result))
    }
    rdata
}
.stop_for_status <-function(response, op) {
        status <- status_code(response)
        if (status < 400L)
            return(invisible(response))

        cond <- http_condition(status, "error")
        type <- headers(response)[["content-type"]]
        msg <- NULL
        if (nzchar(type) && grepl("application/json", type)) {
            content <- as.list(response)
            msg <- content[["message"]]
            if (is.null(msg))
                ## e.g., from bond DRS server
                msg <- content$response$text
        } else if (nzchar(type) && grepl("text/html", type)) {
            ## these pages can be too long for a standard 'stop()' message
            cat(as.character(response), file = stderr())
        }

        message <- paste0(
            "'", op, "' failed:\n  ",
            conditionMessage(cond),
            if (!is.null(msg)) "\n  ", msg
        )
        stop(message, call.=FALSE)
    }
.checkResponse <- function(response) {
    msgs <- response[["messages"]]
    if (!is.null(msgs)) {
        if (grepl("Resource not found", msgs))
            stop(msgs)
        else
            message(response[["messages"]])
    }
    if (!is.null(response[["failedIds"]]))
        warning(
            "IDs not mapped: ",
            paste(response[["failedIds"]], collapse = ", "),
            call. = FALSE
        )
    is.null(response[["results"]])
}
.getResponse <- function(jobId) {
    url <- paste0(.UNIPROT_REST_URL, "idmapping/status/", jobId)
    resp <- GET(url = url)
    content(resp, as = "parsed")
}
.prepQuery <- function(columns, format = "tsv", paginate, pageSize) {
    qlist <- list(format = format)
    if (length(columns))
        qlist <- c(qlist, fields = paste(columns, collapse = ","))
    if (paginate)
        qlist <- c(qlist, size = pageSize)
    qlist
}
.dotter <- function(ndots, maxlength) {
    paste0(
        paste0(rep(".", times = ndots), collapse = ""),
        paste0(rep(" ", times = maxlength-ndots), collapse = ""),
        collapse = ""
    )
}
mapUniProt <- function(
        from = "UniProtKB_AC-ID", to = "UniRef90",
        columns = character(0L), query, verbose = FALSE, debug = FALSE,
        paginate = TRUE, pageSize = 500L
) {
    # stopifnot(
    #     isScalarCharacter(from), isScalarCharacter(to),
    #     isCharacter(query) || is.list(query), isTRUEorFALSE(verbose)
    # )
    if (is.character(query))
        query <- list(ids = paste(query, collapse = ","))
    else if (is.list(query))
        query[["ids"]] <- paste(query[["ids"]], collapse = ",")
    files <- c(query, list(from = from, to = to))
    resp <- httpcache::POST(
        url = .messageDEBUG(paste0(.UNIPROT_REST_URL, "idmapping/run"), debug),
        body = files,
        encode = "multipart",
        accept_json()
    )
    .stop_for_status(resp, "idmapping_run")
    submission <- content(resp, as = "parsed")
    jobId <- submission[["jobId"]]
    if (verbose)
        message("ID Mapping jobId: ", jobId)
    pb <- progress::progress_bar$new(
        format = "  (:spin) waiting for query completion:dots :elapsedfull",
        total = NA, clear = FALSE
    )

    while (.checkResponse(.getResponse(jobId))) {
        for (ndot in seq(0, 10)) {
            pb$tick(tokens = list(dots = .dotter(ndot, 10)))
            Sys.sleep(2/8)
        }
        cat("\n")
    }

    url <- paste0(.UNIPROT_REST_URL, "idmapping/details/", jobId)
    resp <- GET(url = .messageDEBUG(url, debug), accept_json())
    .stop_for_status(resp, "idmapping_details_query")
    details <- content(resp, as = "parsed")
    resurl <- .getResultsURL(details[["redirectURL"]], paginate, debug)
    results <- GET(
        url = resurl,
        query = .prepQuery(columns, pageSize = pageSize, paginate = paginate),
        accept_json()
    )
    .stop_for_status(results, "redirectURL_query")
    .handleResults(results, debug)
}




if (F) {

# Read a gff file
gf <- as_tibble(import(paste0(folder_genomics, "pangenome/", set_name, "/postpanaroo_gffs/", "g4_panaroo.gff"), format = "gff"))
table(gf$eC_number != "None")

gf %>%
    filter(str_detect(name, "hdfR"))





# #top_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name, "/top_gene_fst.csv")) %>% filter(!str_detect(gene, "group"))
# list_tg_fst <- top_gene_fst$gene %>% str_split("~~~") %>% unlist()
# #cat(paste(gene_names, collapse = " "))
# top_gene_or <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv")) %>% filter(!str_detect(gene, "group"))
# list_tg_or <- top_gene_or$gene %>% str_split("~~~") %>% unlist()




if (F) {

# Download NIH COG list
# https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab.txt

cog24 <- read_delim( "~/Downloads/cog-24.def.tab.txt", delim = "\t", col_names = c("COG ID", "COG functional category", "COG name", "Gene name", "Functional pathway", "PubMed ID", "PDB ID"))
cog24 <- cog24 %>% clean_names()

#
set_name = "urbn_mel"
#set_name = "elev_med"
top_gene_fst <- read_csv(paste0(folder_data, "genomics_analysis/fst/", set_name, "/top_gene_fst.csv")) %>% filter(!str_detect(gene, "group"))
list_tg_fst <- top_gene_fst$gene %>% str_split("~~~") %>% unlist()
#cat(paste(gene_names, collapse = " "))
top_gene_or <- read_csv(paste0(folder_data, "genomics_analysis/gene_content/", set_name, "/top_gene_or.csv")) %>% filter(!str_detect(gene, "group"))
list_tg_or <- top_gene_or$gene %>% str_split("~~~") %>% unlist()


# Need the prokka output tsv table for COG id
#genome_id <- paste0("g", c(4,5,6,8,9,11,13,16,17,19))
genome_id <- paste0("g", c(21:27, 31:37, 39, 41, 43))
list_gff <- list()
for (gi in genome_id) {
    list_gff[[gi]] <- read_delim(paste0(folder_data, "genomics/annotation/", gi,"/annotated.tsv"), delim = "\t", skip = 0) %>%
        clean_names() %>%
        rename(cog_id = cog)
}

gffs <- bind_rows(list_gff, .id = "genome_id") #%>% distinct(gene, .keep_all = T)
gffs %>%
    filter(str_detect(gene, "abo")) %>%
    view

cog_tg_fst <- gffs %>%
    filter(gene %in% list_tg_fst) %>%
    left_join(cog24)

table(cog_tg_fst$cog_functional_category)

cog_tg_or <- gffs %>%
    filter(gene %in% list_tg_or) %>%
    left_join(cog24)

table(cog_tg_or$cog_functional_category)










}
}
















}
