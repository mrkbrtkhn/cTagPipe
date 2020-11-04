combinedGuideAnalysis <-function(sumReport,v) {
  print("hallola")
  print(head(sumReport))
  print("running function combinedGuideAnalysis ...")
  dt<-as.data.frame(t(combn(rownames(sumReport),2)))
  dt$sumScores <- apply(dt,1,function(x) return(v$data[[x[1]]]$score+v$data[[x[2]]]$score))
  dt$minScore <- apply(dt,1,function(x) return(min(v$data[[x[1]]]$score,v$data[[x[2]]]$score)))

  dt$sumMutations<-apply(dt,1,function(x) return(sum(sumReport[x[1],]$noOfMutations,sumReport[x[2],]$noOfMutations)))
  dt$combMincoreMMs<-apply(dt,1,function(x) return(min(as.numeric(as.character(sumReport[x[1],]$mincoreMMs)),as.numeric(as.character(sumReport[x[2],]$mincoreMMs)))))
  dt$combMintotalMMs<-apply(dt,1,function(x) return(min(as.numeric(as.character(sumReport[x[1],]$mintotalMMs)),as.numeric(as.character(sumReport[x[2],]$mintotalMMs)))))
  dt$combMincoreMMsExonic<-apply(dt,1,function(x) return(min(as.numeric(as.character(sumReport[x[1],]$mincoreMMsExonic)),as.numeric(as.character(sumReport[x[2],]$mincoreMMsExonic)))))
  dt$combMintotalMMsExonic<-apply(dt,1,function(x) return(min(as.numeric(as.character(sumReport[x[1],]$mintotalMMsExonic)),as.numeric(as.character(sumReport[x[2],]$mintotalMMsExonic)))))

  dt$distBetweenGuides<- apply(dt,1,function(x) return(distance(x=v$data[[x[1]]]$PamRanges,y=v$data[[x[2]]]$PamRanges,ignore.strand=TRUE)))
  dt$guidesOverlap<-apply(dt,1,function(x) return(countOverlaps(v$data[[x[1]]]$guideRanges,  v$data[[x[2]]]$guideRanges,ignore.strand=T)==1))
  dt$sumMutations<-apply(dt,1,function(x) return(sum(sumReport[x[1],]$noOfMutations,sumReport[x[2],]$noOfMutations)))

  dt[order(dt$minScore,decreasing=T),]->dt
  head(dt)
  return(dt)
}
