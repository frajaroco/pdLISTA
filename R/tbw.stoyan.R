#' @title Stoyan rule for temporal bandwidth selection
tbw.stoyan <- function(xyt,t.region,stoyan=0.15){

  if (missing(t.region)){
    t.region <- range(xyt[ ,3],na.rm=TRUE)}
  
  xytimes <- xyt[ ,3]
  npt <- length(xytimes)
  tlength <- max(t.region)-min(t.region)
  
  lambda <- npt/tlength
  # Stoyan & Stoyan 1995, eq (15.16), page 285
  h <- stoyan/sqrt(lambda)
  # conversion to standard deviation
  tbw <- h#/sqrt(5)
  
  invisible(return(tbw)) 
}

