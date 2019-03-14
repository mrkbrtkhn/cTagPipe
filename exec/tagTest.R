setwd("/Users/marek/JLUbox/Andreas/cTagPipe")
library(roxygen2)
library(devtools)
list.files()
install("cTagPipeTest")#,reload = TRUE)
library(cTagPipeTest)
#getMappings
#getMappings("hg38")
#getMappings2("hg38")
#remove.packages("cTagPipeTest")
#detach("package:cTagPipeTest", unload=TRUE)

#Sys.info()["sysname"]->mySys
#if (mySys == "Linux") {source('/data1/COMPUTING/R_data/chipseq_functions.R')}
#if (mySys == "Darwin") {source('/Users/marek/Computing/R_stuff/functions/chipseq_functions.R')}
#library(CRISPRseek)
#library(RCurl)
#library(DT)
#library(msa)
#library(msaR)
#library(seqinr)

#setwd("/Users/marek/Documents/Projects/andreas/tagging")
#source("Crispr_functions_4.R")


###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### target sequences
mAID <- DNAString(c("AAAGAGAAGAGTGCTTGTCCTAAAGATCCAGCCAAACCTCCGGCCAAGGCACAAGTTGTGGGATGGCCACCGGTGAGATCATACCGGAAGAACGTGATGGTTTCCTGCCAAAAATCAAGCGGTGGCCCGGAGGCGGCGGCGTTCGTGAAGGTATCAATGGACGGAGCACCGTACTTGAGGAAAATCGATTTGAGGATGTATAAA"))
tripleFlag<-DNAString(c("GACTACAAGGATGACGACGATAAGGACTACAAGGATGACGACGATAAGGACTACAAGGATGACGACGATAAG"))
t2A<-DNAString(c("GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGACCT"))
fiveP<-DNAString(c("GCATCGTACGCGTACGTGTTTGG"))
threeP<-DNAString(c("CCAAACACGTACGCGTACGATGCG"))
neoM<-DNAString(c("GGATCGGCCATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGGCGATGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCTGA"))
hygM<-DNAString("ATGAAAAAGCCTGAACTCACCGCGACGTCTGTCGAGAAGTTTCTGATCGAAAAGTTCGACAGCGTCTCCGACCTGATGCAGCTCTCGGAGGGCGAAGAATCTCGTGCTTTCAGCTTCGATGTAGGAGGGCGTGGATATGTCCTGCGGGTAAATAGCTGCGCCGATGGTTTCTACAAAGATCGTTATGTTTATCGGCACTTTGCATCGGCCGCGCTCCCGATTCCGGAAGTGCTTGACATTGGGGAATTCAGCGAGAGCCTGACCTATTGCATCTCCCGCCGTGCACAGGGTGTCACGTTGCAAGACCTGCCTGAAACCGAACTGCCCGCTGTTCTGCAGCCGGTCGCGGAGGCCATGGATGCGATCGCTGCGGCCGATCTTAGCCAGACGAGCGGGTTCGGCCCATTCGGACCGCAAGGAATCGGTCAATACACTACATGGCGTGATTTCATATGCGCGATTGCTGATCCCCATGTGTATCACTGGCAAACTGTGATGGACGACACCGTCAGTGCGTCCGTCGCGCAGGCTCTCGATGAGCTGATGCTTTGGGCCGAGGACTGCCCCGAAGTCCGGCACCTCGTGCACGCGGATTTCGGCTCCAACAATGTCCTGACGGACAATGGCCGCATAACAGCGGTCATTGACTGGAGCGAGGCGATGTTCGGGGATTCCCAATACGAGGTCGCCAACATCTTCTTCTGGAGGCCGTGGTTGGCTTGTATGGAGCAGCAGACGCGCTACTTCGAGCGGAGGCATCCGGAGCTTGCAGGATCGCCGCGGCTCCGGGCGTATATGCTCCGCATTGGTCTTGACCAACTCTATCAGAGCTTGGTTGACGGCAATTTCGATGATGCAGCTTGGGCGCAGGGTCGATGCGACGCAATCGTCCGATCCGGAGCCGGGACTGTCGGGCGTACACAAATCGCCCGCAGAAGCGCGGCCGTCTGGACCGATGGCTGTGTAGAAGTACTCGCCGATAGTGGAAACCGACGCCCCAGCACTCGTCCGGAGGCAAAGGAATTCGGGAGATGGGGGAGGCTAACTGAAACACGGAAGGAGACAATACCGGAAGGAACCCGCGCTATGACGGCAATAAAAAGACAGAATAAAACGCACGGGTGTTGGGTCGTTTGTTCATAA")

### user input
### window around stop
winSize<-100
browserWindow<-500
distanceBetweenGuides <- 100
recArmLength<-150
### genome
myGenome<-"hg38"

### elimination rules
# additional perfect hits in same gene
removeAddGeneHits <- TRUE

### set DB
loadTxdb2(myGenome)->txdb
getMappings(txdb)->mGenIDMap
head(mGenIDMap)
#myGene <- "NM_025115"
myGene <- "NM_002106"
myGene<-"NM_052927" # PWWP2A
myGene<-"NM_153252"# 'BRWD3'
#myGene <- c("NM_001304504") # HMG20A


mapRefSeqToSymbol(myGene,mGenIDMap)

### iterate over larger list
mGenIDMap->rfsq
rfsq[grep("NM_",rfsq$tx),]->rfsq
#rfsq[!duplicated(rfsq$gene),]->rfsq
dim(rfsq)
subset(rfsq,gene=="MIER1")
myGene<-"NM_001077700" # MIER1
myGene<-c("NM_000237") ### split stop gene
myGene<-c("NM_001243526")  ### - strand split stop gene
myGene<-c("NM_001243028") ## AKT2 big list test
myGene<-c("NM_001077700") ## MIER1

rfsq[sample(1:nrow(rfsq),100),]->myRfsq
for (i in 1:nrow(myRfsq)) {
  strsplit(as.character(myRfsq$tx[i]),"\\.")[[1]][1]->myGene
  ### get gene symbol
  mapRefSeqToSymbol(myGene,mGenIDMap)->mySymbol
  mySymbol
  print(paste(i,mySymbol))

  extractSeqAroundStop(txdb,myGene,winSize,myGenome)->genRanges
  isPlus<-as.character(strand(genRanges$cds))[1]=="+"

  createSubGtf(genRanges,myGenome,browserWindow)->winTx

  if (!file.exists(paste(myGene,".xls",sep=""))) {
    CCTop(name = myGene, radQ = "single", sequence = as.character(genRanges$sequence), pamType = "NGG", targetLength = "20", sgRNA5 = "NN",
          sgRNA3 = "NN", inVitroTx = "SP6", totalMismatches = "4", useCore = "on", coreLength = "12", coreMismatches = "2",
          species = myGenome, downloadDir = ".")
  } else {print("file exists!")}

  parseCCTopData(paste(myGene,".xls",sep=""))->v
  correctCCTopGeneAnno(v,genRanges,mySymbol) -> v

  firstPassFilterCCResults(v)->v

  evaluateUpCCTop(v,mySymbol,genRanges,winSize,mGenIDMap)->v

  ###

  do.call(rbind,lapply(v$data,bindCCTop))->bindSumReport

  do.call(rbind,lapply(v$data,function(x) return(x[-7])))->sumReport
  createSubGtf(genRanges,myGenome,browserWindow)->winTx
  plotBrowserView(v,genRanges,winTx,browserWindow=browserWindow)
  #plotBrowserView(v,genRanges,winTx=txdb,browserWindow=browserWindow)

  collectMutations(genRanges,v,myGenome)->v

  ###
  summarizeReport(v)->sumReport
  ### clean up
  sumReport<-sumReport[as.logical(sumReport$hasHitinTarget) &
                         (as.logical(unlist(sumReport$pamMutated)) | as.logical(unlist(sumReport$mutated)) | as.logical(unlist(sumReport$pamOverStop)))
                       #& ifelse(!is.na(sumReport$mintotalMMsExonic),as.numeric(as.character(sumReport$mintotalMMsExonic)),TRUE)
                       & !is.na(sumReport$mintotalMMsExonic)
                       ,]

  ### look at data
  combinedGuideAnalysis(sumReport,v)->gDat

  gSelect <- gDat[!gDat$guidesOverlap,][1,]

  extractArmIntervals(gSelect,v,genRanges,myGenome)->armIntervals

  ### mutate homology arms
  mutateSeqs(gSelect,v,genRanges,armIntervals)->armIntervals

  ### combine sequences to rec. template
  if (as.character(strand(genRanges$cds))[[1]]=="+") {
    recTempA <- xscat(fiveP,armIntervals$seqs[[1]],mAID,tripleFlag,t2A,neoM,armIntervals$seqs[[2]],threeP )
    recTempB <- xscat(fiveP,armIntervals$seqs[[1]],mAID,tripleFlag,t2A,hygM,armIntervals$seqs[[2]],threeP )} else {
      recTempA <- xscat(fiveP,reverseComplement(armIntervals$seqs[[2]]),mAID,tripleFlag,t2A,neoM,reverseComplement(armIntervals$seqs[[1]]),threeP )
      recTempB <- xscat(fiveP,reverseComplement(armIntervals$seqs[[2]]),mAID,tripleFlag,t2A,hygM,reverseComplement(armIntervals$seqs[[1]]),threeP )
    }
  write.fasta(recTempA,"test.fa","aa")

  ### check by multiple alignments
  DNAStringSet(c(v$data[[as.character(gSelect[1,1])]]$seq,v$data[[as.character(gSelect[1,2])]]$seq))->guideSeqs
  guideSeqs

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
  guideSeqs

  if (isPlus) {
    myDNASS <- c(DNAStringSet(recTempA),DNAStringSet(genRanges$sequence[1:103]),DNAStringSet(genRanges$sequence[104:203]),
                 DNAStringSet(mAID),DNAStringSet(guideSeqs[[1]]),DNAStringSet(guideSeqs[[2]]),DNAStringSet(fiveP),DNAStringSet(threeP),DNAStringSet(neoM),
                 DNAStringSet(t2A),DNAStringSet(tripleFlag))
  } else {
    myDNASS <- c(DNAStringSet(recTempA),reverseComplement(DNAStringSet(genRanges$sequence[101:203])),reverseComplement(DNAStringSet(genRanges$sequence[1:100])),
                 DNAStringSet(mAID),DNAStringSet(guideSeqs[[1]]),DNAStringSet(guideSeqs[[2]]),DNAStringSet(fiveP),DNAStringSet(threeP),DNAStringSet(neoM),
                 DNAStringSet(t2A),DNAStringSet(tripleFlag))
  }

  ms<-c(1,2,3,4,5,6,7,8,9,10,11)
  ms<-c(1,2,3,5,6)
  myDNASS[ms]->myDNASSms
  names(myDNASSms)<-c("temp","window_with","window_after","mAID","guide_1","guide_2","fiveP","threeP","NeoM","t2A","tFlag")[ms]
  msaClustalOmega(myDNASSms)->seqfile
  msaConvert(seqfile,type="seqinr::alignment")->seqAlign

  write.fasta(as.list(seqAlign$seq),names=seqAlign$nam , file.out="xx.aln")
  msaR("xx.aln", menu=F, overviewbox = F)
  plotBrowserView(v,genRanges,winTx,gSelect,browserWindow=browserWindow)
  
}

##
