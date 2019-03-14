extractArmIntervals <- function(gSelect,v,genRanges,myGenome) {
  print("running function extractArmIntervals ...")
  
  isPlus<-as.character(strand(genRanges$cds))[1]=="+"
  lastbeforeStop <- ifelse(isPlus,start(genRanges$stopCodonPos)-1,start(genRanges$stopCodonPos)-1)
  cutSites<-extractCutSites(gSelect,v)
  rightExtension<-max(cutSites)-(lastbeforeStop+3)
  leftExtension<-lastbeforeStop-min(cutSites)
  
  rightExtension
  leftExtension
  
  if (rightExtension>0) {
    rightLimit<-lastbeforeStop+rightExtension+recArmLength+4
  } else {rightLimit<-lastbeforeStop+recArmLength+4}
  
  if (leftExtension>0) {
    leftLimit<-lastbeforeStop-leftExtension-recArmLength
  } else {    leftLimit<-lastbeforeStop-recArmLength}
  
  intervalRanges<-GRanges(seqnames=seqnames(genRanges$stopWindow),
                          ranges=IRanges(start=c(leftLimit,lastbeforeStop+4),
                                         end=c(lastbeforeStop,rightLimit)),
                          #strand=ifelse(c(isPlus,isPlus),c("+","+"),c("-","-")),
                          #names=ifelse(c(isPlus,isPlus),c("leftArm","rightArm"),c("rightArm","leftArm")))
                          names=c(isPlus,isPlus))
  if (myGenome =="hg38") {
    getSeq(BSgenome.Hsapiens.UCSC.hg38,intervalRanges)->mySeqs
  } else {stop("Genome not known!")}
  
  leftDat<-data.frame(coord=start(intervalRanges)[1]:end(intervalRanges)[1],seq=unlist(strsplit(as.character(mySeqs[[1]]),"")))
  rightDat<-data.frame(coord=start(intervalRanges)[2]:end(intervalRanges)[2],seq=unlist(strsplit(as.character(mySeqs[[2]]),"")))
  return(list(ranges=intervalRanges,seqs=mySeqs,leftDat=leftDat,rightDat=rightDat))
}
