#' @title Product density LISTA functions
#' @description Computes an edge-corrected kernel estimator of the product density LISTA functions.
#' @param xyt Spatial-temporal coordinates \eqn{(x,y,t)} of the point pattern.
#' @param s.region A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
#' @param t.region vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.
#' @param ds A vector of distances \code{u} at which \eqn{\rho^{(2)i}(u,v)} is computed.
#' @param dt A vector of distances \code{v} at which \eqn{\rho^{(2)i}(u,v)} is computed.
#' @param ks A kernel function for the spatial distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} kernel, or \code{"biweight"}.
#' @param kt A kernel function for the temporal distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} kernel, or \code{"biweight"}.
#' @param hs A bandwidth of the kernel function \code{ks}.
#' @details An individual product density LISA functions \eqn{\rho^{(2)i}(.)} should reveal the extent of the contribution of the event \eqn{u_i} to the global estimator of the second-order product density \eqn{\rho^{(2)}(.)}, and may provide a further description of structure in the data (e.g., determining events with similar local structure through dissimilarity measures of the individual LISA functions), for more details see Cressie and Collins (2001).
#' @return A list containing:
#' \itemize{
#'   \item \code{lisa}: A matrix containing the values of the estimation of \eqn{\widehat{\rho}^{(2)i}(r)}  for each one of \eqn{n} points of the process by rows.
#'   \item \code{ds}: If \code{ds} is missing, a vector of distances \code{u} at which \eqn{\rho^{(2)i}(u)} is computed under the restriction \eqn{0<\epsilon<r}.
#'   \item \code{kernel}: A vector of names and bandwidth of the spatial kernel.
#'   \item \code{s.region}: Parameter passed in argument.
#'   }
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.
#' @references Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.
#' @references Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
LISTAfunct <- function(xyt,s.region,t.region,ds,dt,ks="epanech",hs,kt="box",ht,correction=TRUE){
  
  pts <- xyt[, 1:2]
  xytimes <- xyt[ ,3]
  ptsx <- pts[, 1]
  ptsy <- pts[, 2]
  ptst <- xytimes
  npt <- length(ptsx)
  
  if (missing(s.region)){
    W <- ripras(ptsx,ptsy)
    poly <- W$bdry
    X <- poly[[1]]$x
    Y <- poly[[1]]$y
    s.region <- cbind(X,Y)}
  
  if (missing(t.region)) t.region <- range(xyt[ ,3],na.rm=TRUE)
  
  if (missing(hs)){
    sdm <- dist(pts)
    hs <- dpik(sdm,kernel=ks,range.x=c(min(sdm),max(sdm)))}
  
  if (missing(ht)){
    tdm <- dist(ptst)
    ht <- dpik(tdm,kernel=kt,range.x=c(min(tdm),max(tdm)))}
  
  srspt <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  
  if (missing(ds)){
    rect <- as.rectangle(srspt)
    maxs <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs*1.01, maxs, len=50)}
  
  ds <- sort(ds)
  if(ds[1]==0) {
    ds<- ds[-1]}
  
  nds <- length(ds)
  
  if (missing(dt)){
    maxt <- t.region/4
    dt <- seq(ht*1.01, maxt, len=50)}
  
  dt <- sort(dt)
  if(dt[1]==0) {
    dt<- dt[-1]}

  nt <- length(dt)
  
  bsupt <- max(t.region)
  binft <- min(t.region)
  are <- area.owin(s.region)

  kernel=c(ks=ks,hs=hs,kt=kt,ht=ht,edg=correction)

    if (ks=="box") ks=1 	
	else if (ks=="epanech") ks=2
	else if (ks=="biweight") ks=3

      if (kt=="box") kt=1 	
	else if (kt=="epanech") kt=2
	else if (kt=="biweight") kt=3
 
      if (correction==FALSE) edg=0
      else if (correction==TRUE) edg=1
	  
	  hlista <- array(0, dim = c(nds, ndt, npt))
	  
	  for (i in 1:npt){
  		xi <- ptsx[i]
  		yi <- ptsy[i]
	  	ti <- ptst[i]

     listalis <- NULL
	   hlista[,,i] <- localpdf(ptsx,ptsy,ptst,npt,ds,nds,dt,ndt,bsupt,binft,xi,yi,ti,i,ks,kt,hs,ht,are,edg,wrs)
	   }

	invisible(return(list(hlista=hlista,ds=ds,dt=dt,kernel=kernel)))  
}
 