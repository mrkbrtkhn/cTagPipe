extractSeqAroundStop <- function(txdb,myGene,winSize,myGenome) {

  myFilter <- list(tx_name = c(myGene))
  myExons <- exons(txdb, columns=c("EXONID", "TXNAME"),filter=myFilter)
  myTx <- transcripts(txdb,filter=myFilter)
  myCds<-GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
  myTxd<-GenomicFeatures::exonsBy(txdb, by="gene")
  myCds->allCds
  
  myCds <- sort(myCds[[myGene]])
  
  checkForCompleteStop(myCds)->incompleteStop
  
  ### extract genomic sequence +/- winSize bp
  if (as.character(strand(myTx))=="+") {
    stopCodonPos <- GRanges(seqnames=seqnames(myTx),ranges=IRanges(start=end(myCds)[length(myCds)]-2,end=end(myCds)[length(myCds)]))
    stopWindow<-stopCodonPos
    start(stopWindow)<- start(stopWindow)-winSize
    end(stopWindow)<- end(stopWindow)+winSize
  } else {
    stopCodonPos <- GRanges(seqnames=seqnames(myTx),ranges=IRanges(start=start(myCds)[1],end=start(myCds)[1]+2))
    stopWindow<-stopCodonPos
    start(stopWindow)<- start(stopWindow)-winSize
    end(stopWindow)<- end(stopWindow)+winSize
  }
  myTxd[myTxd %over% stopWindow]->myTxd
  if (myGenome =="hg38") {
    getSeq(BSgenome.Hsapiens.UCSC.hg38,stopWindow)->seqs
  } else {stop("Genome not known!")}
  
  allCds[allCds %over% stopWindow]->otherCds
  otherCds[names(otherCds)!=myGene]
  otherCds[names(otherCds)!=myGene]->otherCds
  
  return(list(sequence=seqs[[1]],stopCodonPos=stopCodonPos,stopWindow=stopWindow,cds=myCds,tx=myTxd,otherCds=otherCds))
  
  
}
