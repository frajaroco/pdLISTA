#' @title \eqn{i}-th product density LISTA function
#' @description Computes \eqn{i}-th edge-corrected kernel estimator of the product density LISTA function.
#' @param xyt Spatial-temporal coordinates \eqn{(x,y,t)} of the point pattern.
#' @param i Point position for which you want to calculate a LISTA function.
#' @param s.region A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
#' @param t.region vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.
#' @param ds A vector of distances \code{u} at which \eqn{\rho^{(2)i}(u,v)} is computed.
#' @param dt A vector of distances \code{v} at which \eqn{\rho^{(2)i}(u,v)} is computed.
#' @param ks A kernel function for the spatial distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} kernel, or \code{"biweight"}.
#' @param kt A kernel function for the temporal distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} kernel, or \code{"biweight"}.
#' @param hs A bandwidth of the kernel function \code{ks}.
#' @param ht A bandwidth of the kernel function \code{kt}.
#' @param correction It is \code{TRUE} by default and the Ripley's isotropic edge-correction weights are computed. If it is \code{FALSE} the estimated is without edge-correction.
#' @details An individual product density LISTA functions \eqn{\rho^{(2)i}(.,.)} is a puntual  contribution coming from \eqn{(u_i,t_i)} to the global estimator of the second-order product density \eqn{\rho^{(2)}(.,.)}, and may provide local description of structure in the spatio-temporal point data (e.g., determining events with similar local structure through dissimilarity measures of the individual LISTA functions), for more details see Cressie and Collins (2001).
#' @return A list containing:
#' \itemize{
#'   \item \code{hlista}: A matrix containing the values of the estimation of \eqn{\widehat{\rho}^{(2)i}(r,t)} for the point \eqn{i} of the process by rows.
#'   \item \code{ds}: Vector of distances \code{u} at which \eqn{\rho^{(2)i}(r,t)} is computed under the restriction \eqn{0<\epsilon<r}.
#'   \item \code{dt}: Vector of distances \code{v} at which \eqn{\rho^{(2)i}(r,t)} is computed under the restriction \eqn{0<\delta<t}.
#'   \item \code{kernel}: A vector of names and bandwidth of the spatial and temporal kernel.
#'   }
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.
#' @references Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.
#' @references Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
#' @examples 
#' ## Not run:
#' #################
#' 
#' # Realisations of the homogeneous spatio-temporal Poisson processes
#' stp <- rpp(100)$xyt
#' # Generated spatio-temporal point pattern
#' plot(stp)
#' 
#' # Randomly selected point and its product density LISTA 
#' fixed <-sample(1:length(stp[,1]),1)
#'  
#' # Estimation of the individual i-th product density LISTA function 
#' out <- iLISTA(stp,i=fixed)
#' z1 <- out$hlista
#' 
#' # Spatio-temporal LISTA surface
#' par(mfrow=c(1,1))
#' persp(out$ds,out$dt,z1,theta=-45,phi=30,zlim=range(z1,na.rm=TRUE),expand=0.7,ticktype="detailed",xlab="r = distance",ylab="t = time",zlab="",cex.axis=0.7, cex.lab=0.7)
#' contour(out$ds,out$dt,z1,drawlabels=TRUE,axes=TRUE,xlab="r = distance",ylab="t = time",cex.axis=0.7, cex.lab=0.7)
#' 
#' ## End(Not run)
#' #################
iLISTA <- function(xyt,i,s.region,t.region,ds,dt,ks="epanech",hs,kt="epanech",ht,correction=TRUE){
  
  pts <- xyt[, 1:2]
  xytimes <- xyt[ ,3]
  ptsx <- pts[, 1]
  ptsy <- pts[, 2]
  ptst <- xytimes
  xi <- ptsx[i]
  yi <- ptsy[i]
  ti <- ptst[i]
  npt <- length(ptsx)
  
  if (missing(s.region)){
    W <- ripras(ptsx,ptsy)
    poly <- W$bdry
    X <- poly[[1]]$x
    Y <- poly[[1]]$y
    s.region <- cbind(X,Y)}
  srspt <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  
  if (missing(t.region)) t.region <- range(xyt[ ,3],na.rm=TRUE)
  
  if (missing(hs)){
    sdm <- dist(pts)
    hs <- dpik(sdm,kernel=ks,range.x=c(min(sdm),max(sdm)))}
  
  if (missing(ht)){
    tdm <- dist(ptst)
    ht <- dpik(tdm,kernel=kt,range.x=c(min(tdm),max(tdm)))}
  
  if (missing(ds)){
    rect <- as.rectangle(srspt)
    maxs <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs*1.01, maxs, len=50)}
  
  ds <- sort(ds)
  if(ds[1]==0) {
    ds<- ds[-1]}
  
  nds <- length(ds)
  
  if (missing(dt)){
    maxt <- diff(t.region)/4
    dt <- seq(ht*1.01, maxt, len=50)}
  
  dt <- sort(dt)
  if(dt[1]==0) {
    dt<- dt[-1]}
  
  ndt <- length(dt)
  
  bsupt <- max(t.region)
  binft <- min(t.region)
  are <- area(srspt)
  
  kernel=c(ks=ks,hs=hs,kt=kt,ht=ht,edg=correction)
  
  stpxy <- ppp(x=ptsx,y=ptsy,window=srspt)
  wrs <- edge.Ripley(stpxy,pairdist(pts))
  
  if (ks=="box") ks=1 	
  else if (ks=="epanech") ks=2
  else if (ks=="biweight") ks=3
  
  if (kt=="box") kt=1 	
  else if (kt=="epanech") kt=2
  else if (kt=="biweight") kt=3
  
  if (correction==FALSE) edg=0
  else if (correction==TRUE) edg=1
  
  stlistapd <- array(0, dim = c(nds,ndt)) 
  
  storage.mode(stlistapd) <- "double"
  
  stlistapd <- .Fortran("listapdfunction",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt), 
                        as.double(ds),as.integer(nds),as.double(dt),as.integer(ndt),as.double(bsupt),as.double(binft),
                        as.double(xi),as.double(yi),as.double(ti),as.integer(i),as.integer(ks),as.integer(kt),as.double(hs),
                        as.double(ht),as.double(are),as.integer(edg),as.double(wrs),(stlistapd))
  
  hlista <- stlistapd[[22]]
  
  invisible(return(list(hlista=hlista,ds=ds,dt=dt,kernel=kernel)))  
}