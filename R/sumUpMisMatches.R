sumUpMisMatches<-function(x,mySymbol,genRanges,mGenIDMap) {
  ###
  ### summarize MMs in core and total
  ### works on list of Mismatch data objects
  ###
  for (i in 1:length(x)) {
    x[[i]]$data->dt
    dt$coreMMs <- 0
    dt$totalMMs <- 0
    dt
    for (j in 1:nrow(dt)) {
      strsplit(dt$alignment[j],"\\[")[[1]]->mSplits
      for (k in 1:2) {
        dt$coreMMs[j]<-countPattern("-",mSplits[2])
        dt$totalMMs[j]<-countPattern("-",mSplits[1])+dt$coreMMs[j]
      }
    }
    x[[i]]$data<-dt
    #x[[i]]$hasHitinTarget <- any(grepl(mySymbol,dt$`gene name`))
    x[[i]]$hasHitinTarget <- any(grepl(mySymbol,names(genRanges$tx)))
    #x[[i]]$hitsOtherExonInLocus <-  any(names(genRanges$tx) != mySymbol)
    x[[i]]$hitOverlapsOtherCds <- any(!sapply(names(genRanges$otherCds), function(x) mapRefSeqToSymbol(x,mGenIDMap)==mySymbol))
    x[[i]]$hasHitinOtherExon <- any(grepl("Exonic",dt[dt$`gene name` != mySymbol,]$position))
      
    x[[i]]$mincoreMMs <- min(dt[dt$`gene name` != mySymbol,]$coreMMs)
    x[[i]]$mintotalMMs <-min(dt[dt$`gene name` != mySymbol,]$totalMMs)
    
    if(!any(dt$`gene name` != mySymbol & grepl("Exonic",dt$position))){
      x[[i]]$mincoreMMsExonic <- Inf 
    } else {x[[i]]$mincoreMMsExonic <-  min(dt[dt$`gene name` != mySymbol & grepl("Exonic",dt$position),]$coreMMs)}
   
    if(!any(dt$`gene name` != mySymbol & grepl("Exonic",dt$position))){
      x[[i]]$mintotalMMsExonic <- Inf  
    } else {x[[i]]$mintotalMMsExonic <-  min(dt[dt$`gene name` != mySymbol & grepl("Exonic",dt$position),]$totalMMs)}
    
  }
  return(x)
}
