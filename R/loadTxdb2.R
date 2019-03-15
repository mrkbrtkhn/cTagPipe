loadTxdb2<- function(genome,pathToFile="") {
  if (pathToFile=="") {
    if (genome=="hg38") {
      library(BSgenome.Hsapiens.UCSC.hg38)
      #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
      library(org.Hs.eg.db)

      sysDatFolder <- system.file('extdata',package='cTagPipeTest')
      myGenes<-paste(sysDatFolder,c("hg38_UCSC_genes_iGenomes.gtf"),sep="/")
      myGenesZip<-paste(sysDatFolder,c("hg38_UCSC_genes_iGenomes.gtf.zip"),sep="/")
      gdbFile<-paste(sysDatFolder,c("hg38_UCSC_genes_iGenomes.sqlite"),sep="/")


      if (!file.exists(myGenes) & file.exists(myGenesZip)) {
        unzip(myGenesZip,exdir=sysDatFolder)
      }

      if (!file.exists(gdbFile)) {
        makeTxDbFromGFF(myGenes,format="gtf")->txdb
        saveDb(txdb, file=gdbFile)
      } else {loadDb(gdbFile)->txdb}


    } else {stop("genome not known!")}
  }

  return(txdb)
}
