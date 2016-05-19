localpdf <- function(ptsx,ptsy,ptst,npt,ds,nds,dt,ndt,bsupt,binft,xi,yi,ti,i,ks,kt,hs,ht,are,edg,wrs){

stlistapd <- array(0, dim = c(nds,ndt)) 

storage.mode(stlistapd) <- "double"

stlistapd <- .Fortran("listapdfunction",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt), 
                      as.double(ds),as.integer(nds),as.double(dt),as.integer(ndt),as.double(bsupt),as.double(binft),
                      as.double(xi),as.double(yi),as.double(ti),as.integer(i),as.integer(ks),as.integer(kt),as.double(hs),
                      as.double(ht),as.double(are),as.integer(edg),as.double(wrs),(stlistapd))

hlista <- stlistapd[[22]]

return(hlista=hlista)}