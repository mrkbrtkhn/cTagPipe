checkForCompleteStop <- function(myCds) {
  incompleteStop <- FALSE
  lastExonCodingLength<-NULL
  if (as.character(strand(myCds))[1]=="+") {
    lastExonCodingLength <- width(myCds)[length(myCds)]
    if (lastExonCodingLength < 3) {
      incompleteStop <- TRUE
    }
      
    
  } else {
    lastExonCodingLength <- width(myCds)[1]
    if (lastExonCodingLength < 3) {
      incompleteStop <- TRUE
      }
    }
  return(list(incompleteStop,lastExonCodingLength))
}