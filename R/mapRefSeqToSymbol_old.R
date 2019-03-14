mapRefSeqToSymbol_old <- function(refGene,myMap) {
  
  myMap[[1]]->xx
  myMap[[2]]->yy
  mySymbol<-vector()
  for (i in 1:length(refGene)) {
    xx[refGene[i]]->myEntrez
    yy[[as.character(myEntrez)]]->tmp
    if (i == 1) {tmp->mySymbol} else {append(mySymbol,tmp)->mySymbol}
  }
  return(mySymbol)
}
