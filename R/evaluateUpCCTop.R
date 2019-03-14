evaluateUpCCTop<- function(x,mySymbol,genRanges,winSize,mGenIDMap) {
  print("summarize CCTop results ...")
  sumUpMisMatches(x$data,mySymbol,genRanges,mGenIDMap)->x$data
  print("analyzing guideRNAs overlap with CDS  ...")
  guideOverlapsCodingRegion(x,genRanges,mySymbol,winSize)->x
    
  return(x)
}
