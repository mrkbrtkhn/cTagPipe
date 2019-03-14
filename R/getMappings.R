#' @export
getMappings<-function(txdb){
  as.list(txdb)->txdump
  data.frame(tx=txdump$transcripts$tx_name,gene=txdump$genes$gene_id,stringsAsFactors=FALSE )->rfsq
  #rfsq[grep("NM_",rfsq$tx),]->rfsq
  return(rfsq)
}
