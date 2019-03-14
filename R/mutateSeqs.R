mutateSeqs<-function(gSelect,v,genRanges,armIntervals) {
  print("running function mutateSeqs ...")
 
  
  for (i in c(as.character(gSelect$V1),as.character(gSelect$V2))) {
    #for (i in c("T1","T8")) {
    if (!v$data[[i]]$mutations$pamOverStop) {
      #print(v$data[[i]]$mutations)
      for (j in v$data[[i]]$mutations$mutList) {
        for (k in c(3,4)) {
          if (j$coodinate %in% armIntervals[[k]]$coord) {
            print(paste("mutating position",j$coodinate,"in",names(armIntervals)[k], "to",j$base))
            armIntervals[[k]][armIntervals[[k]]$coord == j$coodinate,]$seq <- j$base
          }
        }
      }
    }
    
    
  }
  
  for (k in c(3,4)) {
    armIntervals$seqs[[k-2]] <- DNAString(paste(armIntervals[[k]]$seq,collapse=""))
  }
  print("done")
  return(armIntervals)
}
