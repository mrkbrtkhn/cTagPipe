summarizeReport<-function(v) {
  print(paste("running summarizeReport function ..."))
  lapply(v$data,function(x) return(x[c(-7,-15,-16,-20)]))->ll
  #str(ll[[1]])
  as.data.frame(do.call(rbind,lapply(ll,unlist)))->sumReport
  as.data.frame(do.call(rbind,lapply(v$data,function(x) return(x$mutations[1:3])))) -> mutReport
  
  unlist(lapply(v$data,function(x) return(length(x$mutations$mutList))))->noOfMutations
  lapply(v$data,function(x) return(lapply(x$mutations$mutList,function(x) return(c(x$type)))))->tmp
  unlist(lapply(tmp,function(x) return(paste(unlist(x),collapse="//"))))->mutType
  cbind(sumReport,mutReport,data.frame(noOfMutations=noOfMutations,mutType=mutType))->dat
  
  return(dat)
}
