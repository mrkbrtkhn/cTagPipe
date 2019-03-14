extractOverlappingExonSeq<-function(genRanges,genome,x) {
  if (genome == "hg38") {require(BSgenome.Hsapiens.UCSC.hg38)
    assign("geno", BSgenome.Hsapiens.UCSC.hg38)} else {stop("genome not known!")}
  cds_seqs <-extractTranscriptSeqs(geno,GRangesList(genRanges$cds[order(genRanges$cds$exon_rank)]) )

  Biostrings::translate(cds_seqs)->aaSeq
  cdsDat<-as.data.frame(genRanges$cds)
  if (nrow(cdsDat)>1) {
  coords<-do.call("c",apply(cdsDat,1,function(x) x[2]:x[3]))
  } else {coords<-cdsDat$start:cdsDat$end}
  if (!any(as.logical(strand(genRanges$cds)=="+"))) {coords<-rev(coords)}
  data.frame(coordinate=coords,
             base=strsplit(as.character(cds_seqs),"")[[1]],
             aa=rep(strsplit(as.character(aaSeq),"")[[1]],each=3),
             aaPos=rep(1:length(strsplit(as.character(aaSeq),"")[[1]]),each=3),
             tripletPos=rep(c(1,2,3),length(strsplit(as.character(aaSeq),"")[[1]])),stringsAsFactors=FALSE)->seqDat
  
  seqDat[1:5,]
  
  ### build splicesite ranges
  startspliceRanges<-GRanges(seqnames=as.character(seqnames(genRanges$stopCodonPos)),ranges=IRanges(start=start(genRanges$cds)-2,end=start(genRanges$cds)-1))
  endspliceRanges<-GRanges(seqnames=as.character(seqnames(genRanges$stopCodonPos)),ranges=IRanges(start=end(genRanges$cds)+1,end=end(genRanges$cds)+2))
  #if (as.character(strand(genRanges$cds))[[1]]=="+") {
    c(startspliceRanges[-1],endspliceRanges[-(length(endspliceRanges))])->spliceRngs
    sort(spliceRngs)->spliceRngs
  #} else {c(startspliceRanges[-(length(startspliceRanges))],endspliceRanges[-1])->spliceRngs}
  
  
  as.data.frame(x$guideRanges)->gDat
  if (gDat$strand == "+") {
    x$seq->gSeq
  } else {reverseComplement(DNAString(x$seq))->gSeq}
  guideDat <- data.frame(coordinate = gDat$start:gDat$end,
                         guideSeq=DNAString(gSeq),
                         guideRevCompSeq=complement (DNAString(gSeq)),stringsAsFactors=FALSE)
  guideDat
  
  merge(seqDat,guideDat, by="coordinate",all=TRUE)->mDat
  mDat[!is.na(mDat$guideSeq),]
  
  
  ### build same structure for PAM sequence
  as.data.frame(x$PamRanges)->pDat
  if (pDat$strand == "+") {
    substr(x$seq,nchar(x$seq)-1,nchar(x$seq))->pSeq
  } else {substr(reverseComplement(DNAString(x$seq)),1,2)->pSeq}
  pamDat <- data.frame(coordinate = pDat$start:pDat$end,
                         pamSeq=DNAString(pSeq),
                         pamRevCompSeq=complement (DNAString(pSeq)),stringsAsFactors=FALSE)
  pamDat
  
  merge(mDat,pamDat, by="coordinate",all=TRUE)->mDat
  mDat[with(mDat,order(mDat$aaPos,mDat$tripletPos)),]->mDat
  #sub select + next AA
  mDat[(min(which(!is.na(mDat$guideSeq)))-2):(max(which(!is.na(mDat$guideSeq)))+2),]->mDat
  
  
  ### build same structure for splice acceptor sites
  as.data.frame(spliceRngs)->ssDat
  ssDat<-data.frame(coordinate=sort(apply(ssDat,1,function(x) return(x[2]:x[3]))),splice=rep("spliceSite",length(sort(apply(ssDat,1,function(x) return(x[2]:x[3]))))))
  merge(mDat,ssDat, by="coordinate",all=TRUE)->mDat
  
  
  ### get genomic sequence
  sSeq<-as.character(getSeq(geno,GRanges(seqnames=as.character(seqnames(genRanges$cds))[1],
                                         ranges=IRanges(start=min(mDat$coordinate[!is.na(mDat$coordinate)]),
                                                        end=max(mDat$coordinate[!is.na(mDat$coordinate)])))))
  
  
  sDat<-data.frame(coordinate=min(mDat$coordinate[!is.na(mDat$coordinate)]):max(mDat$coordinate[!is.na(mDat$coordinate)]),
                   genomeSequence=strsplit(sSeq,"")[[1]],stringsAsFactors=FALSE)
  merge(mDat,sDat, by="coordinate",all=TRUE)->mDat
  if (as.character(strand(genRanges$cds))[[1]]== "+") {
    mDat[order(mDat$coordinate),]->mDat
  } else {
      mDat[order(mDat$coordinate,decreasing=T),]->mDat
           }
  
  
  
  ### clean up NAs
  if (length(mDat[is.na(mDat$aa),]$aa)>0) {
    mDat[is.na(mDat$aa),]$aa<-"none"
    
  }
  
  ### restrict to meaningful window
  mDat[mDat$coordinate>min(mDat$coordinate[!is.na(mDat$guideSeq)])-5 &mDat$coordinate<max(mDat$coordinate[!is.na(mDat$guideSeq)])+5,]->mDat
  
  return(mDat)
}
