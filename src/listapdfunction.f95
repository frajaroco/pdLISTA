!
!     Francisco J. Rodriguez-Cortos, January 2014 
!
!     This function provides an edge-corrected estimate
!     of the spatio-temporal LISTA product density functions.

!

       subroutine listapdfunction(x,y,txy,n,s,ns,t,nt,bsupt,binft, &
      & xi,yi,ti,i,ks,kt,hs,ht,are,edg,wrs,stlistapd)
	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     x,y,txy: coordinates and times of the point process of length n.
!     s: vector of the ns distances at which to calculate the function
!     t: vector of the nt times at which to calculate the function. 
!     bint, bsupt: lower and upper boundaries of the time domain.
!     xi,yi,ti: Fixed coordinate where LISTA function is computed.
!     ks, kt: Spatial and temporal Kernel functons.
!     hs,ht: Spatial and temporal bandwidths.
!     are: Area of the spatial windows.
!     edg: switch on edge-correction factor. 
!     wrs: Ripley's isotropic edge-correction fator.
!     stlistapd: zero matrix of dimension ns x nt.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       implicit real*8(a-h,o-z)

       integer n,ns,nt,iu,iv,ks,kt,i,j,edg
       real*8 stlistapd,two,hs,ht,binft,bsupt,xi,yi,ti,tij,hij,wij,wrs
       real*8 pi,kern,kerns,kernt,are,etij,esij,x,y,xyt,s,t
       dimension x(n),y(n),txy(n),wrs(n,n),s(ns),t(nt),stlistapd(ns,nt)

       two=2d0
       pi=3.14159265d0
	   
            do j=1,n	
            do iu=1,ns
            do iv=1,nt		 
            if (j.ne.i) then
            hij=sqrt(((xi-x(j))*(xi-x(j)))+((yi-y(j))*(yi-y(j))))
            tij=abs(ti-txy(j))
            if (ks.eq.1) then
                kerns=boxkernel((s(iu)-hij)/hs,hs)
                else if (ks.eq.2) then
                    kerns=ekernel((s(iu)-hij)/hs,hs)
                    else if (ks.eq.3) then
                        kerns=qkernel((s(iu)-hij)/hs,hs)
            end if
            if (kt.eq.1) then
                kernt=boxkernel((t(iv)-tij)/ht,ht)
                else if (kt.eq.2) then
                    kernt=ekernel((t(iv)-tij)/ht,ht)
                    else if (kt.eq.3) then
                        kernt=qkernel((t(iv)-tij)/ht,ht)
            end if
            kern=kerns*kernt
            if (kern.ne.0) then
			    if (edg.eq.0) then
                    wij=((n-1)*kern)/(4d0*are*(bsupt-binft)*pi*s(iu))
                    stlistapd(iu,iv)=stlistapd(iu,iv)+wij
                end if
                if (edg.eq.1) then
                    bsup=ti+tij
                    binf=ti-tij
                    if ((bsup.le.bsupt).and.(binf.ge.binft)) then
                        etij=1d0
                    else
		                etij=two
                    end if
             esij=wrs(i,j)
             wij=((n-1)*kern*esij*etij)/(4d0*are*(bsupt-binft)*pi*s(iu))
             stlistapd(iu,iv)=stlistapd(iu,iv)+wij
                end if   
            end if
            end if
            end do
		    end do
            end do 
            
         return

         end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     boxkernel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	     function boxkernel(x,h)

	     implicit real*8 (a-h,o-z)

         double precision x, h

         if (dabs(x).le.1) then
           boxkernel=1d0/2d0
         else
           boxkernel=0d0
         end if
         boxkernel=boxkernel/h

         return
         
         end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Epanechnikov kernel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       function ekernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x

       if (dabs(x).le.1) then
           ekernel=(3d0/4d0)*(1-x**2)
       else
           ekernel=0d0
       end if
       ekernel=ekernel/h

       return
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     quartic (biweight) kernel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
	  function qkernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x, h

       if (dabs(x).le.1) then
           qkernel=(15d0/16d0)*(1-x**2)**2
       else
           qkernel=0d0
       end if
       qkernel=qkernel/h

       return
       end
       
