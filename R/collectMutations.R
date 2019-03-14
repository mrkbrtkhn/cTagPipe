#' Title collectMutations
#'
#' @param genRanges to do
#' @param v results list
#' @param myGenome genome global definition string
#'
#' @return results list v
#' @export
#'
collectMutations<-function(genRanges,v,myGenome) {
  ### iterate through guides and annotates collective mutations
  gList<-list()
  for (i in 1:length(v$data)) {
    x<-v$data[[i]]
    print(i)
    extractOverlappingExonSeq(genRanges,myGenome,x)->exonContext
    selectcMutations(genRanges,exonContext)->v$data[[i]]$mutations
  }
  return(v)
}
