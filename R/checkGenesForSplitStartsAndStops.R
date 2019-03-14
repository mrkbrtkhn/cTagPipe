checkGenesForSplitStartsAndStops<-function(txdb) {
  myCds<-GenomicFeatures::cdsBy(txdb, by="tx", use.names=TRUE)
  startGenes<-list()
  stopGenes<-list()
  for (i in 1:length(myCds)) {
    print(i)
    x<-myCds[[i]]
    width(x)->wth
      if (wth[length(wth)] %in% c(1,2)) {
        stopGenes[[names(myCds)[i]]]<-x
      } 
      if (wth[1] %in% c(1,2)) {
        startGenes[[names(myCds)[i]]]<-x
      }
  }
  return(list(startGenes,stopGenes))
}
