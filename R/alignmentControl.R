alignmentControl <- function(v,gSelect,recTempA,mAID,fiveP,threeP,genRanges,neoM,hygroM,t2A,isPlus) {
  print("running alignmentControl ...")
  
  DNAStringSet(c(v$data[[as.character(gSelect[1,1])]]$seq,v$data[[as.character(gSelect[1,2])]]$seq))->guideSeqs

  
  for (i in 1:2) {
    if (isPlus) {
      if (!as.logical(strand((v$data[[as.character(gSelect[1,i])]]$guideRanges))=="+")) {
        reverseComplement(guideSeqs[[i]])->guideSeqs[[i]]
      }
    } else {
      if (as.logical(strand((v$data[[as.character(gSelect[1,i])]]$guideRanges))=="+")) {
        reverseComplement(guideSeqs[[i]])->guideSeqs[[i]]
      }
    }
  }
  
  print("here we came")
  
  if (isPlus) {
    myDNASS <- c(DNAStringSet(guideSeqs[[1]]),DNAStringSet(recTempA),DNAStringSet(guideSeqs[[2]]),DNAStringSet(genRanges$sequence[1:103]),DNAStringSet(genRanges$sequence[104:203]),
                 DNAStringSet(mAID),DNAStringSet(fiveP),DNAStringSet(threeP),DNAStringSet(neoM),
                 DNAStringSet(t2A),DNAStringSet(tripleFlag))
  } else {
    myDNASS <- c(DNAStringSet(guideSeqs[[1]]),DNAStringSet(recTempA),DNAStringSet(guideSeqs[[2]]),reverseComplement(DNAStringSet(genRanges$sequence[101:203])),reverseComplement(DNAStringSet(genRanges$sequence[1:100])),
                 DNAStringSet(mAID),DNAStringSet(fiveP),DNAStringSet(threeP),DNAStringSet(neoM),
                 DNAStringSet(t2A),DNAStringSet(tripleFlag))
  }
  
  
  
  ms<-c(1,2,3,4,5,6,7,8,9,10,11)
  #ms<-c(1,2)
  myDNASS[ms]->myDNASSms
  names(myDNASSms)<-c("guide_1","rec. template","guide_2","5' genome","3' genome","mAID","fiveP","threeP","NeoM","t2A","tFlag")[ms]
  msaClustalOmega(myDNASSms)->seqfile
  msaConvert(seqfile,type="seqinr::alignment")->seqAlign
#  write.fasta(as.list(seqAlign$seq),names=seqAlign$nam , file.out="xx.aln")
  write.fasta(as.list(seqAlign$seq)[match(names(myDNASSms),seqAlign$nam)],names=names(myDNASSms) , file.out="xx.aln")
  print("done")
}
