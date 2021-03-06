% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iLISA.R
\name{iLISTA}
\alias{iLISTA}
\title{\eqn{i}-th product density LISTA function}
\usage{
iLISTA(xyt, i, s.region, t.region, ds, dt, ks = "epanech", hs,
  kt = "epanech", ht, correction = TRUE)
}
\arguments{
\item{xyt}{Spatial-temporal coordinates \eqn{(x,y,t)} of the point pattern.}

\item{i}{Point position for which you want to calculate a LISTA function.}

\item{s.region}{A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.}

\item{t.region}{vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}

\item{ds}{A vector of distances \code{u} at which \eqn{\rho^{(2)i}(u,v)} is computed.}

\item{dt}{A vector of distances \code{v} at which \eqn{\rho^{(2)i}(u,v)} is computed.}

\item{ks}{A kernel function for the spatial distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} kernel, or \code{"biweight"}.}

\item{hs}{A bandwidth of the kernel function \code{ks}.}

\item{kt}{A kernel function for the temporal distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} kernel, or \code{"biweight"}.}

\item{ht}{A bandwidth of the kernel function \code{kt}.}

\item{correction}{It is \code{TRUE} by default and the Ripley's isotropic edge-correction weights are computed. If it is \code{FALSE} the estimated is without edge-correction.}
}
\value{
A list containing:
\itemize{
  \item \code{hlista}: A matrix containing the values of the estimation of \eqn{\widehat{\rho}^{(2)i}(r,t)} for the point \eqn{i} of the process by rows.
  \item \code{ds}: Vector of distances \code{u} at which \eqn{\rho^{(2)i}(r,t)} is computed under the restriction \eqn{0<\epsilon<r}.
  \item \code{dt}: Vector of distances \code{v} at which \eqn{\rho^{(2)i}(r,t)} is computed under the restriction \eqn{0<\delta<t}.
  \item \code{kernel}: A vector of names and bandwidth of the spatial and temporal kernel.
  }
}
\description{
Computes \eqn{i}-th edge-corrected kernel estimator of the product density LISTA function.
}
\details{
An individual product density LISTA functions \eqn{\rho^{(2)i}(.,.)} is a puntual  contribution coming from \eqn{(u_i,t_i)} to the global estimator of the second-order product density \eqn{\rho^{(2)}(.,.)}, and may provide local description of structure in the spatio-temporal point data (e.g., determining events with similar local structure through dissimilarity measures of the individual LISTA functions), for more details see Cressie and Collins (2001).
}
\examples{
## Not run:
#################

# Realisations of the homogeneous spatio-temporal Poisson processes
stp <- rpp(100)$xyt
# Generated spatio-temporal point pattern
plot(stp)

# Randomly selected point and its product density LISTA 
fixed <-sample(1:length(stp[,1]),1)
 
# Estimation of the individual i-th product density LISTA function 
out <- iLISTA(stp,i=fixed)
z1 <- out$hlista

# Spatio-temporal LISTA surface
par(mfrow=c(1,1))
persp(out$ds,out$dt,z1,theta=-45,phi=30,zlim=range(z1,na.rm=TRUE),expand=0.7,ticktype="detailed",xlab="r = distance",ylab="t = time",zlab="",cex.axis=0.7, cex.lab=0.7)
contour(out$ds,out$dt,z1,drawlabels=TRUE,axes=TRUE,xlab="r = distance",ylab="t = time",cex.axis=0.7, cex.lab=0.7)

## End(Not run)
#################
}
\author{
Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
}
\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
}

