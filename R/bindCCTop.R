bindCCTop <- function(x) {
  return(data.frame(name=rep(x$name,nrow(x$data)),x$data))
}
