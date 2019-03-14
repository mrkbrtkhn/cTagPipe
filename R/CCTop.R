CCTop <- function(name = "test", radQ = "single", sequence = "TTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCA",
                  pamType = "NGG", targetLength = "20", sgRNA5 = "NN",
                  sgRNA3 = "NN", inVitroTx = "T7", totalMismatches = "4", useCore = "on",
                  coreLength = "12", coreMismatches = "2", species = "oryLat2", downloadDir = "."){
  require(httr)
  uri <- "https://crispr.cos.uni-heidelberg.de/cgi-bin/search.py"
  dat.list <- list(name = name,
                   radQ = radQ,
                   sequence = sequence,
                   pamType = pamType,
                   targetLength = targetLength,
                   sgRNA5 = sgRNA5,
                   sgRNA3 = sgRNA3,
                   inVitroTx = inVitroTx,
                   totalMismatches = totalMismatches,
                   useCore = useCore,
                   coreLength = coreLength,
                   coreMismatches = coreMismatches,
                   species = species)

  r <- POST(uri, body = dat.list)
  if(r$status_code == 200){
    id <- unlist(strsplit(r$url, "="))[2] # get id
    downloadLinkResults <- paste0("https://crispr.cos.uni-heidelberg.de/result/", id, "/", dat.list$name, ".xls")
    # create download link for result table
    while(!url.exists(downloadLinkResults)){ # check if url exists
      cat("calculation in progress...\n")
      Sys.sleep(60)
    }
    download.file(downloadLinkResults, destfile = file.path(downloadDir, paste0(dat.list$name, ".xls")))
    # download result table
    downloadLinkFasta <- paste0("https://crispr.cos.uni-heidelberg.de/result/", id, "/", dat.list$name, ".fasta")
    # create download link for result table
    download.file(downloadLinkFasta, destfile = file.path(downloadDir, paste0(dat.list$name, ".fasta")))
    # download result table
  }else{
    error(paste0("Server response is not 200. Response is ", r$status_code, ". Something went wrong!"))
  }
}
