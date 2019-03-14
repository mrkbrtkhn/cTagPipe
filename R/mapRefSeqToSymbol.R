mapRefSeqToSymbol <- function(refGene,myMap) {
  
  if (refGene %in% myMap$tx) {
    mySymbol= myMap$gene[myMap$tx==refGene]
  } else {stop("gene not known!")}
 
  return(mySymbol)
}
