parseCCTopData<-function(cctopFile="test.xls") {
  res<-list()
  readLines(cctopFile)->x
  ### scan header
  grep("^T1",x)[1] -> firstGuideIdx
  for (i in 1:(firstGuideIdx-1)) {
    strsplit(x[i],":")[[1]]->sP
    
    if (!length(sP)==0) {
      if (sP[1]=="Species") {res[["Species"]]<-gsub("\t","",sP[2])}
      if (sP[1]=="Input") {res[["Input"]]<-gsub("\t","",sP[2])}
      if (sP[1]=="PAM") {res[["PAM"]]<-gsub("\t","",sP[2])}
      if (sP[1]=="Target site length") {res[["TargetSiteLength"]]<-as.numeric(gsub("\t","",sP[2]))}
      if (sP[1]=="Target site 5' limitation") {res[["TargetSite5limitation"]]<-(gsub("\t","",sP[2]))}
      if (sP[1]=="Target site 3' limitation") {res[["TargetSite3limitation"]]<-(gsub("\t","",sP[2]))}
      if (sP[1]=="Core length") {res[["CoreLength"]]<-as.numeric(gsub("\t","",sP[2]))}
      if (sP[1]=="Core MM") {res[["CoreMM"]]<-as.numeric(gsub("\t","",sP[2]))}
      if (sP[1]=="Total MM") {res[["TotalMM"]]<-as.numeric(gsub("\t","",sP[2]))}
    }
  }
  
  ### parse off-targets lists
  grep("^T",x)[grep("^T",x)>=firstGuideIdx]->idX
  append(idX, length(x)+1)->idX
  idX
  
  resData<-list()
  expData<-FALSE
  for (i in 2:length(idX)) {
    k<-0
    #  for (i in 2:15) {
    tList<-list()
    n<-1
    for (j in idX[i-1]:(idX[i]-1)) {
      if (n %% 1000 == 0) {print(n)}
      n<-n+1
      strsplit(x[j],"\t")[[1]]->sP
      
      if (!expData) {
        
        if (grepl("^T",sP[1])) {
          tList[["name"]]<-sP[1]
          
          tList[["seq"]]<-sP[2]
          tList[["score"]]<-as.numeric(sP[3])
          tList[["CrispRaterScore"]]<-as.numeric(sP[5])
          
        }
        if (grepl("^Oligo fwd",sP[1])) {tList[["OligoFwd"]]<-sP[2]}
        if (grepl("^Oligo rev",sP[1])) {tList[["OligoRev"]]<-sP[2]}
        #Oligo adding fwd Oligo adding rev
        if (grepl("^Oligo adding fwd",sP[1])) {tList[["OligoFwd"]]<-sP[2]}
        if (grepl("^Oligo adding rev",sP[1])) {tList[["OligoRev"]]<-sP[2]}
        
        if (grepl("^Chromosome",sP[1])) {
          #print(tList[["name"]])
          expData<-TRUE
          sP->myColNames
          as.data.frame(matrix(rep(NA,length(sP)),nrow=1))->myDat
          colnames(myDat)<-sP
        } 
      }
      
      
      ### 
      
      if (expData) {
        k<-k+1
        if (k<8) {
          if (grepl("^Chromosome",x[j-1])) {
            myDat[1,] <- sP
          } 
          if (!grepl("^Chromosome",x[j-1]) & !grepl("^Chromosome",x[j])) {rbind(myDat,sP)->myDat}
        }
        
        if (length(sP)==0) {
          #print(sP)
          tList[["data"]]<-myDat
          expData<-FALSE
        }
      }
    }
    print(paste(i,j))
    resData[[tList[["name"]]]]<-tList
  }
  res[["data"]]<-resData
  return(res)
}