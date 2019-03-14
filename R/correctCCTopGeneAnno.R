correctCCTopGeneAnno <- function(v,genRanges,mySymbol) {
  for (i in 1:length(v$data)) {
    v$data[[i]]$data->dt
    gRangT <- GRanges(seqnames=dt$Chromosome,ranges=IRanges(start=as.numeric(dt$start),end=as.numeric(dt$end)))
    dt$'gene name'[ which(gRangT %over% genRanges$stopWindow)] <- mySymbol
  }
  return(v)
}
