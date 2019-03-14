createSubGtf<-function(genRanges,myGenome,browserWindow) {
  Sys.info()["sysname"]->mySys
  if (mySys == "Linux") {
    myBedTools <- paste(system.file('Ubuntu64',package='cTagPipeTest'),"/bedtools",sep="")
  }
  if (mySys == "Darwin") {
    myBedTools <- paste(system.file('OSX',package='cTagPipeTest'),"/bedtools",sep="")
  }
  

  if (myGenome=="hg38") {
    myGenes<-c("hg38_UCSC_genes_iGenomes.gtf")
  }
  
  exRanges<-genRanges$stopWindow
  start(exRanges)<-start(exRanges)-1000
  end(exRanges)<-end(exRanges)+1000

  export.bed(exRanges,con="stopWindow.bed")

  if (myGenome=="hg38") {
    sysDatFolder <- system.file('extdata',package='cTagPipeTest')
    myGenes<-paste(sysDatFolder,c("hg38_UCSC_genes_iGenomes.gtf"),sep="/")
  } 
  mycall <- paste(myBedTools,"intersect -wa -a",myGenes ,"-b stopWindow.bed > window.gtf")
  write.table(mycall,file="subTx.sh",sep="",quote=F,row.names=F,col.names=F)
  system("sh subTx.sh")
#  system("sed 's/stop_codon/CDS/g' window.gtf > tmp")
#  system("mv tmp window.gtf")
  
#  makeTxDbFromGFF("window.gff3",format="gff3")->winTx
  makeTxDbFromGFF("window.gtf",format="gtf")->winTx
  return(winTx)

}
