
extractCutSites<- function(gSelect,v){
  x<-list()
  x[[1]]<-ifelse(as.character(strand(v$data[[as.character(gSelect[,1])]]$PamRanges))=="+",
                 start(v$data[[as.character(gSelect[,1])]]$PamRanges)-3,
                 end(v$data[[as.character(gSelect[,1])]]$PamRanges)+3)
  x[[2]]<-ifelse(as.character(strand(v$data[[as.character(gSelect[,2])]]$PamRanges))=="+",
                 start(v$data[[as.character(gSelect[,2])]]$PamRanges)-3,
                 end(v$data[[as.character(gSelect[,2])]]$PamRanges)+3)
  return(unlist(x))
}
