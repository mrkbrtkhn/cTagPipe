loadTxdb<- function(genome) {

  if (genome=="hg38") {
  
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    
    myGenes<-c("genes/hg38_UCSC_genes_iGenomes.gtf")
    gdbFile<-c("genes/hg38_UCSC_genes_iGenomes.sqlite")
    
    if (!file.exists(gdbFile)) {
      makeTxDbFromGFF(myGenes,format="gtf")->txdb
      saveDb(txdb, file=gdbFile)
    } else {loadDb(gdbFile)->txdb}
  } else {stop("genome not known!")}
  return(txdb)
}
