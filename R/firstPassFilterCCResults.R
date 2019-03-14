 firstPassFilterCCResults<- function(v) {
  elimList<-list()
  for (i in 1:length(v$data)) {
    ### check for additional perfect hits
    if (removeAddGeneHits) {
      subset(v$data[[i]]$data,MM==0)->x
      if (nrow(x)>1) {
        #print(paste("guide",i,""))
        if (duplicated(x$`gene name`)[2]) {
          elimList<-c(elimList,i)
        }
      }
      
    }
  }
  if (length(elimList)>0) {
    v$data <- v$data[-unlist(elimList)]
  }
  return(v)
  
  }
