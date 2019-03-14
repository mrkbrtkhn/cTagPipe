#' @export
getMappings_old<-function(myGenome){
  if (myGenome=="hg38") {
    print("mapping genes to symbols")
    
    library(org.Hs.eg.db)
    ### get ENTREZ ID
    x <- org.Hs.eg.db::org.Hs.egREFSEQ2EG
    ### get corresponding SYMBOL
    y <- org.Hs.eg.db::org.Hs.egSYMBOL
    # Get the RefSeq identifier that are mapped to an entrez gene ID
    mapped_seqs <- AnnotationDbi::mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_seqs])
    
    
    
    
    # Get the entrez gene identifiers that are mapped to a gene symbol
    mapped_genes <- AnnotationDbi::mappedkeys(y)
    # Convert to a list
    yy <- as.list(y[mapped_genes])
  }
  return(list(xx,yy))
}