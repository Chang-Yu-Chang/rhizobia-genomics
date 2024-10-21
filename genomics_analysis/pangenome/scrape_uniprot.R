#' Scrape UnitProt to convert gene names into accession IDs

library(httr)

isJobReady <- function(jobId) {
    pollingInterval = 5
    nTries = 20
    for (i in 1:nTries) {
        url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
        r <- GET(url = url)
        status <- content(r, as = "parsed", encoding = "UTF-8")
        if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
            return(TRUE)
        }
        if (!is.null(status[["messages"]])) {
            print(status[["messages"]])
            return (FALSE)
        }
        Sys.sleep(pollingInterval)
    }
    return(FALSE)
}

getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
        url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
        url <- gsub("/results/", "/results/stream/", redirectURL)
    }
}

files = list(
    taxon_id = "266834",
    ids = "fixN,nifH",
    from = "Gene_Name",
    to = "UniProtKB"
)
files = list(
    taxon_id = "266834",
    ids = "P21802,P12345",
    from = "UniProtKB_AC-ID",
    to = "UniRef90"
)

r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())
submission <- content(r, as = "parsed")
submission

if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable = read.table(text = content(r), sep = "\t", header=TRUE)
    print(resultsTable)
}

