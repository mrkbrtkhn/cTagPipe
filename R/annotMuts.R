annotMuts<-function(genRanges,myGenome,x) {
  extractOverlappingExonSeq(genRanges,myGenome,x)->exonContext
  exonContext$mutations<-NA
  if (!x$mutations$pamOverStop) {
    for (i in x$mutations$mutList) {
      exonContext[which(exonContext$coordinate==as.numeric(i[[2]])),]$mutations<-i[[3]]
    } 
  }
  
  return(exonContext)
}
