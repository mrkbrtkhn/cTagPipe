plotBrowserView<-function(v,genRanges,winTx,gSelect=NULL,browserWindow) {
  makeGRangesFromDataFrame(do.call(rbind,lapply(v$data,function(x) return(x$data[1,c(1:5)]))))->gR
  gR$guide<-names(gR)
  gR
  mycol=rep(8,length(v$data))
  if (!is.null(gSelect)) {
    mycol[which(names(v$data) %in% c(as.character(gSelect$V1),as.character(gSelect$V2)))]<-c(3,4)
    print(mycol)
  }

  grTrack<-GeneRegionTrack(winTx,chromosome = as.character(seqnames(genRanges$stopWindow)),
                           start=start(genRanges$stopWindow)-browserWindow,
                           end=end(genRanges$stopWindow)+browserWindow,
                           showId = TRUE,geneSymbols=TRUE,name = "Gene Model",
                           col.sampleNames="black",background.title="white",col.axis="black",
                           col.title="black",col=1,fill=1,
                           fontfamily="Helvetica")
  annoTrack1 <-AnnotationTrack(gR,background.title="white",col.title="black",name="guideRNAs",
                               featureAnnotation = "guide",fill=mycol)
  axisTrack <- GenomeAxisTrack()
  plotTracks(list(grTrack,axisTrack,annoTrack1),chromosome = as.character(seqnames(genRanges$stopWindow)),
             from=start(genRanges$stopWindow)-browserWindow,
             to=end(genRanges$stopWindow)+browserWindow)
}

