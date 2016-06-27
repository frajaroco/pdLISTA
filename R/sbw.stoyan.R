#' @title Stoyan rule for spatial bandwidth selection
sbw.stoyan <- function(xyt,s.region,stoyan=0.15){

  pts <- xyt[,1:2]
  ptsx <- pts[, 1]
  ptsy <- pts[, 2]
  npt <- length(pts[,1])
  
if (missing(s.region)){
  W <- ripras(ptsx,ptsy)
  poly <- W$bdry
  X <- poly[[1]]$x
  Y <- poly[[1]]$y
  s.region <- cbind(X,Y)}
  srspt <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  
  areaW <- area(srspt)
  lambda <- npt/areaW
  
  # Stoyan & Stoyan 1995, eq (15.16), page 285
  h <- stoyan/sqrt(lambda)
  # conversion to standard deviation
  sbw <- h#/sqrt(5)
  
  invisible(return(sbw=sbw)) 
}
