guideOverlapsCodingRegion<-function(xfull,genRanges,mySymbol,winSize) {
  for (i in 1:length(xfull$data)) {
    xfull$data[[i]]->x
    #print(paste("processing",x[[i]]$name))
    matchPattern(x$seq,genRanges$sequence)->mfx
    #print(x$seq)
    #print(mfx)
    if (length(nchar(mfx))>0) {
      mfx -> mx
      mStrand<-"+"
    }
    matchPattern(reverseComplement(DNAString(x$seq)) ,genRanges$sequence)->mrx
    if (length(nchar(mrx))>0) {
      mrx -> mx
      mStrand<-"-"}
    start(mx)
    end(mx)
    ###guide coordinates
    start(genRanges$stopCodonPos)-winSize+start(mx)-1->gStart
    start(genRanges$stopCodonPos)-winSize+end(mx)-1->gEnd
    ###PAM coordinates
    ifelse(mStrand=="+",gEnd-1,gStart)->PamStart
    ifelse(mStrand=="+",gEnd,gStart+1)->PamEnd

    x$guideRanges<-GRanges(seqnames=seqnames(genRanges$stopCodonPos),
                           ranges=IRanges(start=gStart,
                                          end=gEnd), strand=mStrand,name=x$name)
    x$PamRanges <-GRanges(seqnames=seqnames(genRanges$stopCodonPos),
                          ranges=IRanges(start=PamStart,
                                         end=PamEnd), strand=mStrand,name=x$name)
    x$guideOverCds<-overlapsAny(x$guideRanges ,genRanges$cds,ignore.strand=TRUE) 
    x$PamOverCds<-overlapsAny(x$PamRanges, genRanges$cds,ignore.strand=TRUE)
    
    ### report distance to stop
    x$distancePamToStop<-distance(x$PamRanges,genRanges$stopCodonPos)

    
    x->xfull$data[[i]]
  }
  
  return(xfull)
}
